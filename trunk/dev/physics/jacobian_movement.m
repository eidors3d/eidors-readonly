function J = jacobian_movement(fwd_model, img)
% JACOBIAN_MOVEMENT   Computes the Jacobian matrix for conductivity and
% electrode movement variables in 3D EIT.
% Args:     fwd_model - the EIDORS object forward model
%            img - the image background conductivity
%
% fwd_model.conductivity_jacobian - function to calculate conductivity
%                                   Jacobian (defaults to jacobian_adjoint)
%
% Returns:          J - the Jacobian matrix [Jc, Jm]
%
% WARNING: THIS CODE IS EXPERIMENTAL AND GIVES PROBLEMS
% SEE: Camille Gomez-Laberge, Andy Adler
% Direct EIT Jacobian calculations for conductivity change
%  and electrode movement,  Physiol. Meas., 29:S89-S99, 208

% (C) 2007, Camille Gomez-Laberge and Andy Adler.
%  License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return ; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_MOVEMENT with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
   warning off EIDORS:DeprecatedInterface

end
fwd_model= img.fwd_model;

img = physics_data_mapper(img);
if ~ismember(img.current_physics, supported_physics)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'JACOBIAN_MOVEMENT',img.current_physics);
end

% System matrix and its parameters
pp = fwd_model_parameters( fwd_model );
pp.dfact = factorial(pp.n_dims);
pp.DEBUG = 0;
if pp.DEBUG
    pp.ss_mat = calc_unconnected_system_mat( fwd_model, img);
    pp.fwd_meas =fwd_solve( fwd_model, img);
end

%Tensor product for conductivity matrix %MC 25/05/2012
I_nd=speye(pp.n_dims+1); 
sigma_mat=spdiags(img.elem_data,0,pp.n_elem,pp.n_elem);
pp.kron_cond=kron(sigma_mat,I_nd);

pp.Ce= connectivity_matrix( pp );
s_mat= calc_system_mat( img );
[pp.Vc, pp.Re] = Vc_Re_matrices( pp, fwd_model, s_mat.E );

ws = warning('query','EIDORS:OverwritingData');
warning off EIDORS:OverwritingData
if isfield(fwd_model,'conductivity_jacobian')
   img.fwd_model.jacobian = fwd_model.conductivity_jacobian;
   Jc= calc_jacobian( img );
else
   img.fwd_model = mdl_normalize(fwd_model, 0); % we normalize on our own
   Jc = jacobian_adjoint(img);
%    Jc = calc_conductivity_jacobian(pp, fwd_model, img);
end
warning(ws.state,'EIDORS:OverwritingData');


Jm = calc_movement_jacobian(pp, fwd_model, img);
J=[Jc,Jm];

if pp.normalize
    data= fwd_solve( img );
    J= J ./ (data.meas(:)*ones(1,size(J,2)));
end



% Define the element connectivity matrix Ce
function Ce= connectivity_matrix( pp );
lengthX = (pp.n_dims+1)*pp.n_elem;
%%%%%lengthY = pp.n_node; %MC 11/05/2012
lengthY=size(pp.N2E,2);
Xidx = pp.ELEM(:);
Yidx = ones(lengthX, 1);
Ce = sparse(1:lengthX, Xidx, Yidx, lengthX, lengthY);



% Calculate fwd_solution and Impedance mapper matrices
function [Vc, Re] = Vc_Re_matrices( pp, fwd_model, s_mat );
% Define the stimulation matrix Vc for each node I, and pattern
% Ground node is never excited; remove it from the index
%%%%%%nodeidx = 1:pp.n_node; %MC 11/05/2012
nodeidx = 1:size(s_mat);
nodeidx( fwd_model.gnd_node ) = [];

% The stimulation matrix Vc is the voltage at each node (row) for a
% stimulation (column)
Vc = zeros(pp.n_node, pp.n_stim);
Vc(nodeidx, :) = s_mat(nodeidx, nodeidx) \ pp.QQ(nodeidx,:);

% Define the electrode resistance matrix Re
% Calculate the resistance between electrodes (row) and all nodes (col)
% N2E matrix maps each electrode to its node(s); we exclude GND
Re = zeros(pp.n_elec, pp.n_node);
Re(:, nodeidx) = pp.N2E(:, nodeidx) / s_mat(nodeidx, nodeidx);

% FIXME: why do we calculate the negative??
Re = -Re;



% Calculate Meas jacobian from derivative on nodes
% Input delVc
% Ouput J
function J= nodes_to_stim_jacobian( delV, fwd_model, pp )

sz_out= size(delV,3);
% Define the conductivity Jacobian Jc
J = zeros(pp.n_meas, sz_out);
% Calculate the Jacobian columns for each stimulation pattern
rowidx = 0;
for j = 1:pp.n_stim
    % Get the measurement pattern for the stimulation pattern j
    meas_pattern = fwd_model.stimulation(j).meas_pattern;
    n_measj = size(meas_pattern, 1);
    % Extract the voltage sensitivity for electrode j
    delVcj = reshape( delV(:,j,:), pp.n_elec, sz_out);
    % Calculate sensitivity block for measurements during stimulation j
    J(rowidx+(1:n_measj), :) = meas_pattern*delVcj;
    rowidx = rowidx+n_measj;
end



% CONDUCTIVITY JACOBIAN (Based on Andy Adler's 1996 algorithms)
% Define the voltage sensitivity delVc on electrode I, for stimulation
% pattern J, for a change in conductivity of element K as a 3D array
function Jc = calc_conductivity_jacobian(pp, fwd_model, img);
delVc = zeros(pp.n_elec, pp.n_stim, pp.n_elem);
% Calculate the values for the voltage sensitivity for each element
for k = 1:pp.n_elem
    if ~mod(k,500)
        fprintf('   JC: element # %d\n',k);
    end
    % Extract the coordinates of the element's four nodes
    Ae = pp.NODE(:,pp.ELEM(:,k))';
    % Augment Ae by adding a column of ones to invert
    Ae = inv([ones(pp.n_dims+1,1), Ae]);
    % Define Be as the matrix Ae with row 1 deleted
    Be = Ae(2:pp.n_dims+1,:);
    % Calculate the system submatrix subSe for the element i
    %%%%%subSe = 2*Be'*Be/pp.dfact/abs(det(Ae)); %MC 11/05/2012
    subSe = Be'*Be/pp.dfact/abs(det(Ae)); 
    % Select the same submatrix of Ce
    subidx = (pp.n_dims+1)*(k-1)+1 : (pp.n_dims+1)*k;
    % The system submatrix is given by the product
    if ~pp.DEBUG
        delVc(:,:,k) = pp.Re * pp.Ce(subidx,:)' * subSe * pp.Ce(subidx,:)...
            * pp.Vc;
    else
        sz= (pp.n_dims+1)*pp.n_elem;
        delSe = sparse(sz,sz);
        se_idx= (pp.n_dims+1)*k+(-pp.n_dims : 0);
        delSe(se_idx, se_idx) = subSe;

        delVc(:,:,k) = pp.Re * pp.Ce' * delSe * pp.Ce * pp.Vc;

        if mod(k,50) == 0
            delta=1e-6;
            img_delta = img;
            img_delta.elem_data(k) = img_delta.elem_data(k) + delta;
            ss_mat_delta= calc_unconnected_system_mat( fwd_model, img_delta);
            delSe_pert = (ss_mat_delta - pp.ss_mat) / delta;

            if norm(delSe -delSe_pert ,1) > 1e-6
                eidors_msg('delSe calc wrong',1);
            end
        end
    end

end
Jc= nodes_to_stim_jacobian( delVc, fwd_model, pp );



% MOVEMENT JACOBIAN
function Jm = calc_movement_jacobian(pp, fwd_model, img)
% The movement Jacobian is defined for each coordinate Jm = [Jmx Jmy Jmz]
% Define the voltage sensitivity delVm on electrode I, for stimulation
% pattern J, for a movement of electrode K as a 3D array
delVm = zeros(pp.n_elec, pp.n_stim, pp.n_elec*pp.n_dims);

% Precalculate
Re_Ce      = pp.Re * pp.Ce';
cond_Ce_Vc = pp.kron_cond * pp.Ce * pp.Vc;
for colidx = 1:pp.n_dims
    fprintf('   JM: direction # %d\n',colidx);
    % Calculate the values for the voltage sensitivity for each electrode
    for k = 1:pp.n_elec
        % Find which elements touch this electrode
        elec_nodes = fwd_model.electrode(k).nodes;
 
        %for compound electrodes, average jacobian for each node
%         delVm_part = zeros(pp.n_elec,pp.n_stim);
        delVm_part = calc_delVm(elec_nodes(:)',pp,fwd_model,img,colidx,...
                Re_Ce, cond_Ce_Vc);
%         delVm_part = sparse(length(pp.kron_cond),length(pp.kron_cond));
%         for each_elec_node= elec_nodes(:)';
%            delVm_part =  delVm_part + ...
%                 calc_delVm(each_elec_node,pp,fwd_model,img,colidx,...
%                 Re_Ce, cond_Ce_Vc);
%         end
%         delVm_part = Re_Ce * delVm_part * cond_Ce_Vc;
        %%%%% delVm_part = delVm_part/length(elec_nodes); %MC 25/05/2012
%         delVm_part = delVm_part;

        vm_idx= k + pp.n_elec*(colidx-1);
        delVm(:,:,vm_idx) = delVm_part;

        if pp.DEBUG
            delta=1e-8;
            mdl_delta = fwd_model;
            mdl_delta.nodes(elec_nodes, colidx) = ...
                mdl_delta.nodes(elec_nodes, colidx) + delta;
            img_delta = img;
            img_delta.fwd_model = mdl_delta;
            S= calc_system_mat(img_delta); 
            S=S.E;
% FIXME: AA+CG 30/1/12
            [Vc_delta] = Vc_Re_matrices( pp, mdl_delta, S);
            delVm_pert = pp.N2E*(Vc_delta - pp.Vc) / delta;
            nn = norm(delVm_part - delVm_pert,1 ); % WHY NEGATIVE?

            %%%%% if nn > 5e-5 ; keyboard; end %MC 25/05/2012
            if nn > 5e-3 ; keyboard; end

        end
    end
end
Jm= nodes_to_stim_jacobian( delVm, fwd_model, pp );



function delVm=  calc_delVm( elec_nodes_array, pp, fwd_model, img, colidx, Re_Ce, cond_Ce_Vc)
I = []; J=[]; S= [];
for elec_nodes= elec_nodes_array(:)';
[rowidx, elemidx] = find(pp.ELEM == elec_nodes);
% Define the system sensitivity matrix to movement delSm
sz= (pp.n_dims+1)*pp.n_elem;
%delSm = sparse(sz,sz);
% For each touching element, calculate the perturbation
jcount = 1;
for j = elemidx'
    % Extract the coordinates of the element's four nodes
    Ae = pp.NODE(:,pp.ELEM(:, j))';
    % Define the invertible matrix P: augment Ae by adding a
    % column of ones to invert
    P = [ones(pp.n_dims+1,1), Ae];
    Ae = inv(P);
    absdetAe = abs(det(Ae));
    % Define Be as the matrix Ae with row 1 deleted
    Be = Ae(2:pp.n_dims+1,:);
    % For this coordinate, perturb P by [rowidx,colidx], which are
    % our paper's perturbation vectors [a,b]
    a = zeros(pp.n_dims+1,1);
    b = a;
    a(rowidx(jcount)) = 1;
    jcount = jcount + 1;
    b(colidx+1) = 1;
    % Calculate the system submatrix subSm for the element j by
    % asymmetric perturbation of the electrode node k
    deldetAe =   1/absdetAe*b'*Ae*a;
    delBe = -Ae*a*b'*Ae;
    delBe = delBe(2:pp.n_dims+1,:);
    %%%%%subSm = 2/pp.dfact*(...
    %%%%%    deldetAe*Be'*Be + ...
    %%%%%    delBe'*Be/absdetAe + ...
    %%%%%    Be'*delBe/absdetAe); %MC 11/05/2012
    subSm = 1/pp.dfact*(...
        deldetAe*Be'*Be + ...
        delBe'*Be/absdetAe + ...
        Be'*delBe/absdetAe);

    if pp.DEBUG
        delta=1e-8;
        %%%%% subSe = 2*Be'*Be/pp.dfact/abs(det(Ae)); %MC 25/05/2012
        subSe = Be'*Be/pp.dfact/abs(det(Ae));
        d_NODE= pp.NODE;
        d_NODE(colidx,elec_nodes) =  d_NODE(colidx,elec_nodes) + delta;
        Ae = d_NODE(:,pp.ELEM(:, j))';
        Ae = inv( [ones(pp.n_dims+1,1), Ae] );
        absdetAe_pert = abs(det(Ae));
        deldetAe_pert = (absdetAe_pert - absdetAe) / delta;
        % Define Be as the matrix Ae with row 1 deleted
        Be = Ae(2:pp.n_dims+1,:);
        %%%%% subSe_delta = 2*Be'*Be/pp.dfact/abs(det(Ae)); %MC 25/05/2012
        subSe_delta = Be'*Be/pp.dfact/abs(det(Ae));
        subSm_pert= (subSe_delta - subSe ) / delta;
        if norm(subSm_pert - subSm,1) > 1e-5
            eidors_msg('subSm calc wrong',1);
            dd= (subSm_pert - 2/pp.dfact/absdetAe *...
                (delBe'*Be + Be'*delBe) )./(Be'*Be);

            fprintf('colidx=%d, j=%d std=%6.4f >',...
                colidx,j, std(dd(:)));
            keyboard
            subSm= subSm_pert;
        end
    end

    % Embed subSm into delSm such that subSm(1,1) is the
    % (4j+1,4j+1) element of delSm
    se_idx= (pp.n_dims+1)*j+(-pp.n_dims : 0);
    switch pp.n_dims
       case 2
          Iidx = vertcat(se_idx,se_idx,se_idx);
          I = [I Iidx(:)];
          J = [J, se_idx,se_idx,se_idx];
          
       case 3
          Iidx = vertcat(se_idx,se_idx,se_idx,se_idx);
          I = [I Iidx(:)];
          J = [J, se_idx,se_idx,se_idx,se_idx];
    end
    S = [S subSm(:)];
%     delSm(se_idx, se_idx) = subSm;
end
end
delSm = sparse(I,J,S,sz,sz);
delVm = Re_Ce * delSm * cond_Ce_Vc;


% The system submatrix is given by the product where delSm is
% non-zero only in submatrices corresponding to touching elements
%%%%% delVm = pp.Re * pp.Ce' * delSm * pp.Ce * pp.Vc; %MC 25/05/2012
% delVm = Re_Ce * delSm * cond_Ce_Vc;
if pp.DEBUG
    delta=1e-8;
    mdl_delta = fwd_model;
    mdl_delta.nodes(elec_nodes, colidx) = ...
        mdl_delta.nodes(elec_nodes, colidx) + delta;
    ss_mat_delta= calc_unconnected_system_mat( mdl_delta, img );
     delSm_pert = (ss_mat_delta - pp.ss_mat) / delta;
    % delSe_pert shound be Ce'*delSe*Ce
    if norm(delSm -delSm_pert ,1) > 1e-5
        eidors_msg('delSm calc wrong',1);
        %%%%% delVm = pp.Re * pp.Ce' * delSm_pert * pp.Ce * pp.Vc; %MC 25/05/2012
        delVm = pp.Re * pp.Ce' * delSm_pert * pp.kron_cond * pp.Ce * pp.Vc;
    
        keyboard
    end
end



function SS= calc_unconnected_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

p= fwd_model_parameters( fwd_model );

d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

SSiidx= floor([0:d*e-1]'/d)*d*ones(1,d) + ones(d*e,1)*(1:d) ;
SSjidx= [1:d*e]'*ones(1,d);
SSdata= zeros(d*e,d);
dfact= (d-1)*(d-2); % Valid for d<=3
for j=1:e
    a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
    idx= d*(j-1)+1 : d*j;
    SSdata(idx,1:d)= a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));

end %for j=1:ELEMs
idx= 1:e*d;
SS= sparse(SSiidx,SSjidx,SSdata) * ...
    sparse(idx,idx, img.elem_data(ceil(idx/d)) );


function do_unit_test;
   unit_test_compare_approaches
   unit_test_matrix_derivatives
   unit_test_diff_jacobian_b2C_const_cond
   unit_test_diff_jacobian_n3r2_const_cond
   unit_test_diff_jacobian_b2C_rand_cond
   unit_test_diff_jacobian_n3r2_rand_cond
  unit_test_3d_inv_solve1
  unit_test_3d_inv_solve2
   
function unit_test_compare_approaches
   inv_model = mk_common_model('d2t2',16);
   img  = mk_image( inv_model);
   fwd_model = inv_model.fwd_model;

   pp = fwd_model_parameters( fwd_model );
   pp.DEBUG = 0;
   pp.dfact = factorial(pp.n_dims);
   s_mat= calc_system_mat(img );
   [pp.Vc, pp.Re] = Vc_Re_matrices( pp, fwd_model, s_mat.E );
   pp.Ce= connectivity_matrix( pp );

   Jc1= calc_conductivity_jacobian(pp, fwd_model, img);
   Jc2= jacobian_adjoint(fwd_model,img);
   unit_test_cmp('Compare J d2t2', Jc1, Jc2, 1e-13);

   fwd_model.normalize_measurements = 1;
   unit_test_cmp('Compare J norm', Jc1, Jc2, 1e-13);



% TEST CODE FOR MATRIX DERIVATIVES
function unit_test_matrix_derivatives

TEST= 'd/dt det(X + t*a*b'')';
d= 1e-8;
X= rand(5);
a=zeros(5,1); a(ceil(5*rand))=1;
b=zeros(5,1); b(ceil(5*rand))=1;
dX_p= (det(X + d*a*b') - det(X) )/d;
dX  = b'*inv(X)*a*det(X);
unit_test_cmp(TEST, dX_p,dX,1e-6);

TEST= 'd/dt inv(X + t*a*b'')';
dX_p= (inv(X + d*a*b') - inv(X) )/d;
dX  = -inv(X)*a*b'*inv(X);
unit_test_cmp(TEST, dX_p,dX,1e-5);

% TEST d/dt 1/abs(det(X + t*a*b')) = abs(1/det(X+t*a*b'))
for i=1:10
    TEST = sprintf('d/dt abs(1/det(X+t*a*b'')) [%02d]',i);
    X= rand(5);
    a=zeros(5,1); a(ceil(5*rand))=1;
    b=zeros(5,1); b(ceil(5*rand))=1;
    dX_p= (1/abs(det(X + d*a*b')) - 1/abs(det(X)) )/d;
    dX  = - 1/abs(det(X))*b'*inv(X)*a;
    unit_test_cmp(TEST, norm([dX_p-dX]),0, 1e-5*norm(dX));
end


   
function unit_test_diff_jacobian_b2C_const_cond
   TEST= 'J_perturb-J_direct - b2C model (const sigma)';
   mdl3dim = mk_common_model( 'b2C' );
   img = mk_image(mdl3dim);
   J_pert=jacobian_movement_perturb(img);
   J_direct =jacobian_movement(img);
   unit_test_cmp(TEST, norm([J_pert-J_direct]),0, 1e-5*norm(J_direct));

   
function unit_test_diff_jacobian_n3r2_const_cond
   TEST= 'J_perturb-J_direct - n3r2 model (const sigma)';
   mdl3dim = mk_common_model( 'n3r2', [16,2] );
   img = mk_image(mdl3dim);
   J_pert=jacobian_movement_perturb(img);
   J_direct =jacobian_movement(img);
   unit_test_cmp(TEST, norm([J_pert-J_direct]),0, 1e-5*norm(J_direct));
   
function unit_test_diff_jacobian_b2C_rand_cond
   TEST= 'J_perturb-J_direct - b2C model (rand sigma)';
   mdl3dim = mk_common_model( 'b2C' );
   cond=0.5+rand(size(mdl3dim.fwd_model.elems,1),1);
   img = mk_image(mdl3dim,cond);
   J_direct =jacobian_movement(img);
   J_pert=jacobian_movement_perturb(img);   
   unit_test_cmp(TEST, norm([J_pert-J_direct]),0, 1e-5*norm(J_direct));
   
function unit_test_diff_jacobian_n3r2_rand_cond
   TEST= 'J_perturb-J_direct - n3r2 model (rand sigma)';
   mdl3dim = mk_common_model( 'n3r2', [16,2] );
   cond=0.5+rand(size(mdl3dim.fwd_model.elems,1),1);
   img = mk_image(mdl3dim,cond);
   J_pert=jacobian_movement_perturb(img);
   J_direct =jacobian_movement(img);
   unit_test_cmp(TEST, norm([J_pert-J_direct]),0, 1e-5*norm(J_direct));
  
function unit_test_3d_inv_solve1
   mdl3dim = mk_common_model( 'n3r2', [16,2] );
   img = mk_image(mdl3dim);
   vh = fwd_solve( img );
   mdl3dim.fwd_model.jacobian = @jacobian_movement;

   mdl3dim.RtR_prior = @prior_movement;

   imgM = inv_solve(mdl3dim, vh, vh);

 function unit_test_3d_inv_solve2
    fmdl = mk_library_model('adult_male_16el_lungs');
    [fmdl.stimulation, fmdl.meas_select] = ...
       mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
    
    outline = shape_library('get','adult_male','boundary');
    
    minnode = min(fmdl.nodes);
    maxnode = max(fmdl.nodes);
    imgsz = [32 32];
    
    xgrid = linspace(minnode(1),maxnode(1),imgsz(1));
    ygrid = linspace(minnode(2),maxnode(2),imgsz(2));
    rmdl = mk_grid_model([],xgrid,ygrid);
    
    % remove pixels outside the model
    x_avg = conv2(xgrid, [1,1]/2,'valid');
    y_avg = conv2(ygrid, [1,1]/2,'valid');
    [x,y] = ndgrid( x_avg, y_avg);
    inside = inpolygon(x(:),y(:),outline(:,1),outline(:,2));
    ff = find(~inside);
    
    rmdl.elems([2*ff, 2*ff-1],:)= [];
    rmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
    rmdl.coarse2fine(:,ff)= [];
    rmdl.mk_coarse_fine_mapping.f2c_offset = [0 0 0.5];
    rmdl.mk_coarse_fine_mapping.z_depth = 0.25;
    
    
    % calculate coarse2fine
    fmdl.coarse2fine = mk_coarse_fine_mapping(fmdl,rmdl);
    
    imdl = select_imdl( fmdl,{'Basic GN dif'});
    imdl.rec_model = rmdl;
    
    img = mk_image(imdl,1);
   vh = fwd_solve( img );
   imdl.prior_use_fwd_not_rec = 1;
   imdl.fwd_model.jacobian = @jacobian_movement;
   imdl.RtR_prior = @prior_movement;
   imgM = inv_solve(imdl, vh, vh);
 
   
   function str = supported_physics
    str = {'conductivity'};
