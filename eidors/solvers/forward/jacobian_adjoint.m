function J= jacobian_adjoint( fwd_model, img)
% JACOBIAN_ADJOINT: J= jacobian_adjoint( img ) 
% Calculate Jacobian Matrix for current stimulation EIT
% J         = Jacobian matrix
% img.fwd_model = forward model
% img.elem_data = background for jacobian calculations
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements


% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_ADJOINT with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end

img = data_mapper(img);
if ~ismember(img.current_params, supported_params)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'JACOBIAN_ADJOINT',img.current_params);
end

org_params = img.current_params;
% all calcs use conductivity
img = convert_img_units(img, 'conductivity');

img.elem_data = check_elem_data(img);

fwd_model= img.fwd_model;

pp= fwd_model_parameters( fwd_model, 'skip_VOLUME' );
s_mat= calc_system_mat( img );

d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

idx= 1:size(s_mat.E,1);
idx( fwd_model.gnd_node ) = [];

sv= zeros(n, pp.n_stim );
sv( idx,:) = left_divide(s_mat.E(idx,idx) , pp.QQ( idx,: ));

zi2E= zeros(pp.n_elec, n);
% the minus below used to be missing
zi2E(:, idx)= -pp.N2E(:,idx)/ s_mat.E(idx,idx) ;

FC= system_mat_fields( fwd_model );

if isfield(fwd_model,'coarse2fine') && strcmp(org_params, 'conductivity');
   DE = jacobian_calc(pp, zi2E, FC, sv, fwd_model.coarse2fine);
   nparam= size(fwd_model.coarse2fine,2);
else
   DE = jacobian_calc(pp, zi2E, FC, sv);
   nparam= e;
end

J = zeros( pp.n_meas, nparam );
idx=0;
for j= 1:pp.n_stim
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   DEj = reshape( DE(:,j,:), pp.n_elec, nparam );
   J( idx+(1:n_meas),: ) = meas_pat*DEj;
   idx= idx+ n_meas;
end

if ~strcmp(org_params,'conductivity')
    J = apply_chain_rule(J, img, org_params);
    if isfield(fwd_model, 'coarse2fine') && ...
          size(img.elem_data,1)==size(fwd_model.coarse2fine,1)
            J=J*fwd_model.coarse2fine;
            nparam = size(fwd_model.coarse2fine,2);
    end
end

%restore img to original condition
if ~strcmp(org_params,'conductivity') || isfield(img, org_params)
    img = rmfield(img,'elem_data');
    img.current_params = [];
end

% calculate normalized Jacobian
if pp.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,nparam));
   
end

% This was a correction for the missing minus above
% J= -J;

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
function DE = jacobian_calc(pp, zi2E, FC, sv, c2f)
d= pp.n_dims+1;
dfact= factorial(d);

do_c2f = ( nargin==5 );
zi2E_FCt = zi2E * FC';

FC_sv   = FC * sv;

if ~do_c2f
   DE= zeros(pp.n_elec, pp.n_stim, pp.n_elem);
   for k= 1:pp.n_elem
       idx= (d-1)*(k-1)+1 : (d-1)*k;
       dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
       DE(:,:,k)= dq;
   end
else
   DE= zeros(pp.n_elec, pp.n_stim, size(c2f,2) );
   if 0 % Code is slower
      de= pp.n_elem * (d-1);
      for k= 1:size(c2f,2);
          chg_col = kron( c2f(:,k), ones(d-1,1));
          dDD_dEj = spdiags(chg_col,0, de, de);
          dq= zi2E_FCt * dDD_dEj * FC_sv;
          DE(:,:,k)= dq;
      end
   else
      de= pp.n_elem * (d-1);
      for k= 1:size(c2f,2);
          ff = find( c2f(:,k) );
          lff= length(ff)*(d-1);
          ff1= ones(d-1,1) * ff(:)';
          ffd= (d-1)*ff1 + (-(d-2):0)'*ones(1,length(ff));
          dDD_dEj = spdiags(c2f(ff1,k), 0, lff, lff);
          dq= zi2E_FCt(:,ffd) * dDD_dEj * FC_sv(ffd,:);
          DE(:,:,k)= dq;
      end
   end
end

function J = apply_chain_rule(J, img, org_params)

switch(org_params)
    case 'resistivity'
        dCond_dPhys = -img.elem_data.^2;
    case 'log_resistivity'
        dCond_dPhys = -img.elem_data;
    case 'log_conductivity'
        dCond_dPhys = img.elem_data;
    otherwise
        error('not implemented yet')
end

J = J.*repmat(dCond_dPhys ,1,size(J,1))';

function elem_data = check_elem_data(img)
   elem_data = img.elem_data; 
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('jacobian_adjoin: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(img.fwd_model, 'coarse2fine');
     c2f = img.fwd_model.coarse2fine;
     sz_c2f = size(c2f);
     switch sz_elem_data(1)
       case sz_c2f(1); % Ok     
       case sz_c2f(2); elem_data = c2f * elem_data;
       otherwise; error(['jacobian_adjoint: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(img.fwd_model)
       error(['jacobian_adjoint: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end



function str = supported_params
    str = {'conductivity'
           'resistivity'
           'log_conductivity'
           'log_resistivity'};
      

function do_unit_test
   current = 4; measure=1;
   [R,img] = test_2d_resistor(current,measure);
   img.fwd_solve.get_all_nodes = 1;
   vs = fwd_solve_1st_order( img);
   va= measure*current*sum(R); % analytic
   unit_test_cmp('2D resistor test', va, vs.meas, 1e-12);

   J = jacobian_adjoint(img);
   unit_test_cmp('2D resistor Jacobian', size(J), ...
      [length(img.fwd_model.stimulation), size(img.fwd_model.coarse2fine,2)]);
   unit_test_cmp('2D resistor Jacobian', std(J),0, 1e-12);
%  unit_test_cmp('2D R voltages', vs.volt(1:3:10)-vs.volt(1), ...

   [R,img] = test_3d_resistor(current,measure);
   img.fwd_solve.get_all_nodes = 1;
   vs = fwd_solve_1st_order( img);
   va= current*sum(R);
   unit_test_cmp('3D resistor test', va, vs.meas, 1e-10);
   J = jacobian_adjoint(img);
   unit_test_cmp('3D resistor Jacobian', size(J), ...
      [length(img.fwd_model.stimulation), size(img.fwd_model.coarse2fine,2)]);
   unit_test_cmp('3D resistor Jacobian', std(J),0, 1e-12);

function [R,img] = test_2d_resistor(current,measure)
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 3; len = 12; 

   fmdl=mk_grid_model([],linspace(0,wid,3), linspace(0,len,4));
   fmdl.electrode(1).nodes = find(fmdl.nodes(:,2) ==   0);
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,2) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R = len / wid / conduc;
   Contact_R = z_contact/wid;
   R = [Block_R, 2*Contact_R];

function [R,img] = test_3d_resistor(current,measure);;
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 2; len = 5; hig=3; 

   fmdl=mk_grid_model([],0:wid, 0:hig, 0:len);
   fmdl.electrode(1).nodes = find(fmdl.nodes(:,3) ==   0);
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,3) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R =  len / wid / hig / conduc;
   Contact_R = z_contact/(wid*hig);
   R = [Block_R, 2*Contact_R];
