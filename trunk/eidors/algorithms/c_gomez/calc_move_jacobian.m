function J = calc_move_jacobian(fwd_model, img_bkgd)
% CG_CALC_JACOBIAN   Computes the Jacobian matrix for conductivity and
% electrode movement variables in 3D EIT.
% Note: this Jacobian is NOT normalized
% Args:     fwd_model - the EIDORS object forward model
%            img_bkgd - the image background conductivity
% Returns:          J - the Jacobian matrix [Jc, Jm]

% (C) 2007, Camille Gomez-Laberge and Andy Adler.
%  License: GPL version 2 or version 3
% $Id$

% System matrix and its parameters

pp = aa_fwd_parameters( fwd_model );
pp.dfact = factorial(pp.n_dims);
pp.DEBUG = 0;
if pp.DEBUG
    pp.ss_mat = calc_unconnected_system_mat( fwd_model, img_bkgd);
    pp.fwd_meas =fwd_solve( fwd_model, img_bkgd);
end

pp.Ce= connectivity_matrix( pp );
s_mat= calc_system_mat( fwd_model, img_bkgd );
[pp.Vc, pp.Re] = Vc_Re_matrices( pp, fwd_model, s_mat.E );
Jc = calc_conductivity_jacobian(pp, fwd_model, img_bkgd);
Jm = calc_movement_jacobian(pp, fwd_model, img_bkgd);
J=[Jc,Jm];

if pp.normalize
    data= fwd_solve( img_bkgd );
    J= J ./ (data.meas(:)*ones(1,size(J,2)));
end



% Define the element connectivity matrix Ce
function Ce= connectivity_matrix( pp );
lengthX = (pp.n_dims+1)*pp.n_elem;
lengthY = pp.n_node;
Xidx = pp.ELEM(:);
Yidx = ones(lengthX, 1);
Ce = sparse(1:lengthX, Xidx, Yidx, lengthX, lengthY);



% Calculate fwd_solution and Impedance mapper matrices
function [Vc, Re] = Vc_Re_matrices( pp, fwd_model, s_mat );
% Define the stimulation matrix Vc for each node I, and pattern
% Ground node is never excited; remove it from the index
nodeidx = 1:pp.n_node;
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
function Jc = calc_conductivity_jacobian(pp, fwd_model, img_bkgd);
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
    subSe = 2*Be'*Be/pp.dfact/abs(det(Ae));
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
            img_delta = img_bkgd;
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
function Jm = calc_movement_jacobian(pp, fwd_model, img_bkgd)
% The movement Jacobian is defined for each coordinate Jm = [Jmx Jmy Jmz]
% Define the voltage sensitivity delVm on electrode I, for stimulation
% pattern J, for a movement of electrode K as a 3D array
delVm = zeros(pp.n_elec, pp.n_stim, pp.n_elec*pp.n_dims);
for colidx = 1:pp.n_dims
    fprintf('   JM: direction # %d\n',colidx);
    % Calculate the values for the voltage sensitivity for each electrode
    for k = 1:pp.n_elec
        % Find which elements touch this electrode
        elec_nodes = fwd_model.electrode(k).nodes;
 
        %for compound electrodes, average jacobian for each node
        delVm_part = zeros(pp.n_elec);
        for each_elec_node= elec_nodes(:)';
           delVm_part =  delVm_part + ...
                calc_delVm(each_elec_node,pp,fwd_model,img_bkgd,colidx);
        end
        delVm_part = delVm_part/length(elec_nodes);

        vm_idx= k + pp.n_elec*(colidx-1);
        delVm(:,:,vm_idx) = delVm_part;

        if pp.DEBUG
            delta=1e-8;
            mdl_delta = fwd_model;
            mdl_delta.nodes(elec_nodes, colidx) = ...
                mdl_delta.nodes(elec_nodes, colidx) + delta;
            [Vc_delta] = Vc_Re_matrices( pp, mdl_delta, ...
                calc_system_mat( mdl_delta, img_bkgd));
            delVm_pert = pp.N2E*(Vc_delta - pp.Vc) / delta;
            nn = norm(delVm_part - delVm_pert,1 ); % WHY NEGATIVE?

            if nn > 5e-5 ; keyboard; end
        end
    end
end
Jm= nodes_to_stim_jacobian( delVm, fwd_model, pp );



function delVm=  calc_delVm( elec_nodes, pp, fwd_model, img_bkgd, colidx)
[rowidx, elemidx] = find(pp.ELEM == elec_nodes);
% Define the system sensitivity matrix to movement delSm
sz= (pp.n_dims+1)*pp.n_elem;
delSm = sparse(sz,sz);
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
    subSm = 2/pp.dfact*(...
        deldetAe*Be'*Be + ...
        delBe'*Be/absdetAe + ...
        Be'*delBe/absdetAe);

    if pp.DEBUG
        delta=1e-8;
        subSe = 2*Be'*Be/pp.dfact/abs(det(Ae));
        d_NODE= pp.NODE;
        d_NODE(colidx,elec_nodes) =  d_NODE(colidx,elec_nodes) + delta;
        Ae = d_NODE(:,pp.ELEM(:, j))';
        Ae = inv( [ones(pp.n_dims+1,1), Ae] );
        absdetAe_pert = abs(det(Ae));
        deldetAe_pert = (absdetAe_pert - absdetAe) / delta;
        % Define Be as the matrix Ae with row 1 deleted
        Be = Ae(2:pp.n_dims+1,:);
        subSe_delta = 2*Be'*Be/pp.dfact/abs(det(Ae));
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
    delSm(se_idx, se_idx) = subSm;
end

% The system submatrix is given by the product where delSm is
% non-zero only in submatrices corresponding to touching elements
delVm = pp.Re * pp.Ce' * delSm * pp.Ce * pp.Vc;
if pp.DEBUG
    delta=1e-8;
    mdl_delta = fwd_model;
    mdl_delta.nodes(elec_nodes, colidx) = ...
        mdl_delta.nodes(elec_nodes, colidx) + delta;
    ss_mat_delta= calc_unconnected_system_mat( mdl_delta, img_bkgd );
    delSm_pert = (ss_mat_delta - pp.ss_mat) / delta;
    % delSe_pert shound be Ce'*delSe*Ce
    if norm(delSm -delSm_pert ,1) > 1e-5
        eidors_msg('delSm calc wrong',1);
        delVm = pp.Re * pp.Ce' * delSm_pert * pp.Ce * pp.Vc;
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

p= aa_fwd_parameters( fwd_model );

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
    SSdata(idx,1:d)= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));
end %for j=1:ELEMs
idx= 1:e*d;
SS= sparse(SSiidx,SSjidx,SSdata) * ...
    sparse(idx,idx, img.elem_data(ceil(idx/d)) );

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST CODE FOR MATRIX DERIVATIVES

% TEST dertiv of det(X + t*a*b')

d= 1e-6;
X= rand(5);
a=zeros(5,1); a(ceil(5*rand))=1;
b=zeros(5,1); b(ceil(5*rand))=1;
dX_p= (det(X + d*a*b') - det(X) )/d;
dX  = b'*inv(X)*a*det(X);
disp([dX_p dX]);

% TEST dertiv of inv(X + t*a*b')
dX_p= (inv(X + d*a*b') - inv(X) )/d;
dX  = -inv(X)*a*b'*inv(X);
disp(norm([dX_p-dX],1));

% TEST dertiv of 1/abs(det(X + t*a*b')) = abs(1/det(X+t*a*b'))
for i=1:20
    X= rand(5);
    a=zeros(5,1); a(ceil(5*rand))=1;
    b=zeros(5,1); b(ceil(5*rand))=1;
    dX_p= (1/abs(det(X + d*a*b')) - 1/abs(det(X)) )/d;
    dX  = - 1/abs(det(X))*b'*inv(X)*a;
    disp(norm([dX_p-dX])/norm(dX));
end
