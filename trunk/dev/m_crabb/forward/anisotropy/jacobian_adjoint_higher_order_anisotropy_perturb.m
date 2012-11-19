function J= jacobian_adjoint_higher_order_anisotropy_perturb( fwd_model, img)
% JACOBIAN_PERTURB: J= jacobian_perturb( fwd_model, img)
% Calculate Jacobian Matrix, based on small perturbations
%   in the forward model. This will tend to be slow, but
%   should be best used to 'sanity check' other code
%
% J         = Jacobian matrix
% fwd_model = forward model
% fwd_model.jacobian_perturb.delta   - delta perturbation to use
% img = image background for jacobian calc

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: jacobian_perturb.m 3326 2012-07-01 18:09:20Z bgrychtol $

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return;end

if isfield(fwd_model,'jacobian_perturb')
   delta = fwd_model.jacobian_perturb.delta;
else
   delta= 1e-6; % tests indicate this is a good value
end

%Get number elements, dimension and CDOF for Jacobian (degrees per element)
n_elem = size(fwd_model.elems,1);
node_dim = size(fwd_model.nodes,2);
CDOF = 0.5*node_dim*(node_dim+1);

% force image to use provided fwd_model
img.fwd_model= fwd_model;

% solve one time to get the size
d0= fwd_solve( img );

%Array for the Jacobian
J= zeros(size(d0.meas,1), n_elem*CDOF);

%2D and 3D are fundamentally different now
if(node_dim==2)
    for i=1:n_elem
        J(:,(i-1)*CDOF +1 )= perturb(img, i, 1, 1, delta, d0);
        J(:,(i-1)*CDOF +2 )= perturb(img, i, 1, 2, delta, d0);
        J(:,(i-1)*CDOF +3 )= perturb(img, i, 2, 2, delta, d0);
    end        
else
    for i=1:n_elem
        J(:,(i-1)*CDOF +1 )= perturb(img, i, 1, 1, delta, d0);
        J(:,(i-1)*CDOF +2 )= perturb(img, i, 1, 2, delta, d0);
        J(:,(i-1)*CDOF +3 )= perturb(img, i, 1, 3, delta, d0);
        J(:,(i-1)*CDOF +4 )= perturb(img, i, 2, 2, delta, d0);
        J(:,(i-1)*CDOF +5 )= perturb(img, i, 2, 3, delta, d0);
        J(:,(i-1)*CDOF +6 )= perturb(img, i, 3, 3, delta, d0);   
    end
end


function Jcol= perturb( img, i, a, b, delta, d0)
   img.elem_data(i,1,a,b)= img.elem_data(i,1,a,b) + delta;
   di= fwd_solve( img );
   Jcol = (1/delta) * (di.meas - d0.meas);

%{
function do_unit_test
% Perturbation Jacobians
% $Id: jacobian_perturb.m 3326 2012-07-01 18:09:20Z bgrychtol $

% imdl= mk_common_model('c2c2',16);
  imdl= mk_common_model('a3cr',16);
% imdl= mk_common_model('n3r2',16);
  imdl.fwd_model.nodes = imdl.fwd_model.nodes*.25;
  img= calc_jacobian_bkgnd(imdl);

  img.fwd_model = mdl_normalize(img.fwd_model, 0);

  img.fwd_model.jacobian=   @jacobian_adjoint;
  img.fwd_model.system_mat= @system_mat_1st_order;
  img.fwd_model.solve=      @fwd_solve_1st_order;
  J_aa= calc_jacobian( img ); % 2 for bug in my code

  img.fwd_model.jacobian=   @jacobian_perturb;
  img.fwd_model.system_mat= @system_mat_1st_order;
  img.fwd_model.solve=      @fwd_solve_1st_order;
  J_aa_p= calc_jacobian( img ); % 2 for bug in my code

  img.fwd_model.jacobian=   @np_calc_jacobian;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np= calc_jacobian( img );

  img.fwd_model.jacobian=   @jacobian_perturb;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np_p= calc_jacobian( img );

  tol = 4e-4;
  unit_test_cmp('J_aa - J_aa_p', norm(J_aa - J_aa_p,'fro')/norm(J_aa,'fro'), 0, tol);
  unit_test_cmp('J_np - J_np_p', norm(J_np - J_np_p,'fro')/norm(J_np,'fro'), 0, tol);
  unit_test_cmp('J_np - J_aa_p', norm(J_np - J_aa_p,'fro')/norm(J_np,'fro'), 0, tol);
  unit_test_cmp('J_np - J_aa  ', norm(J_np - J_aa  ,'fro')/norm(J_np,'fro'), 0, tol);

% Demo model with EIDORS3D
  imdl= mk_common_model('n3r2',16);
  img= calc_jacobian_bkgnd(imdl);
  for i=1:16; imdl.fwd_model.electrode(i).z_contact= .01;end

% Calculate the normal Jacobian
  img.fwd_model.jacobian=   @np_calc_jacobian;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np= calc_jacobian( img );

function display_jacobian_deltas
% Calculate the perturbation Jacobian
  vh= fwd_solve(img);

  for el= 1:50:size(img.fwd_model.elems,1);
     fprintf('\nelem#%03d: ',el);
     for delta= [1e-4,1e-6,1e-8] % delta is perturb
        img_i = img;
        img_i.elem_data(el)= img_i.elem_data(el) + delta;
        vi= fwd_solve(img_i);
        J_p = (vi.meas - vh.meas)/delta; % perturb Jacobian
        fprintf(' %8.6f', norm(J_p - J_np(:,el))/norm(J_np(:,el)));
     end
  end
%}
