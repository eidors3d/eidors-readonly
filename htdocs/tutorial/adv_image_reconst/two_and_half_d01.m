% FIRST RUN DEMO_REAL

if ~exist('demo_img');
   [inhomg_img, demo_img] = demo_real;
end

% Create 2D FEM of all NODES with z=0
f_mdl = demo_img.fwd_model;
n2d = f_mdl.nodes( ...
           (f_mdl.nodes(:,3) == 0), 1:2);
e2d = delaunayn(n2d);
c_mdl = eidors_obj('fwd_model','2d','elems',e2d,'nodes',n2d);
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl );

% Simulate images
vi= fwd_solve(inhomg_img);
homg_img= inhomg_img; homg_img.elem_data(:) = 1;
vh= fwd_solve(homg_img);


imdl.name= 'Nick Polydorides EIT inverse';
imdl.solve=       @np_inv_solve;
imdl.hyperparameter.value = 1e-2;
imdl.R_prior= @np_calc_image_prior;
imdl.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
imdl.jacobian_bkgnd.value= 1;
imdl.reconst_type= 'difference';
imdl.fwd_model= demo_img.fwd_model;


imdl.coarse_fine.solve = imdl.solve;
imdl.solve = @coarse_fine_solve;
% imdl.coarse_fine.mapping = speye(n_e3d); % original solver
  imdl.coarse_fine.mapping = c2f;
imdl= eidors_obj('inv_model', imdl);

img= inv_solve(imdl, vh, vi);
show_slices(img,4);
