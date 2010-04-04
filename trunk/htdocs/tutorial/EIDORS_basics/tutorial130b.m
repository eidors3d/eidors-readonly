% Compare 3D algorithms
% $Id$

% Create Inverse Model
inv3d= eidors_obj('inv_model', 'EIT inverse');
inv3d.reconst_type= 'difference';
inv3d.jacobian_bkgnd.value = 1;
inv3d.fwd_model= imb.fwd_model;
inv3d.fwd_model.np_fwd_solve.perm_sym= '{y}';

% Nick Polydorides' Gauss-Newton Solver
inv3d.hyperparameter.value = 1e-3;
inv3d.solve=       @np_inv_solve;

% Tikhonov prior
inv3d.R_prior=     @tikhonov_image_prior;
imgr(1)= inv_solve( inv3d, vh, vi);
imgn(1)= inv_solve( inv3d, vh, vi_n);

% Nick Polydorides' Prior (Laplace)
inv3d.R_prior=     @np_calc_image_prior;
inv3d.np_calc_image_prior.parameters= [3 1]; %  deg=1, w=1
imgr(2)= inv_solve( inv3d, vh, vi);
imgn(2)= inv_solve( inv3d, vh, vi_n);

% Andrea Borsic's PDIPM TV solver
inv3d.ab_calc_tv_prior.alpha2 = 1e-5;
inv3d.parameters.max_iterations= 20;
inv3d.parameters.term_tolerance= 1e-3;
inv3d.R_prior=     @ab_calc_tv_prior;
inv3d.solve=       @ab_tv_diff_solve;

imgr(3)= inv_solve( inv3d, vh, vi);
imgn(3)= inv_solve( inv3d, vh, vi_n);

% Output image
posn= [inf,inf,2.5,1,1;inf,inf,1.5,1,2;inf,inf,0.5,1,3];
clf;
show_slices(imgr, posn);
print_convert('tutorial130b.png', '-density 100')
show_slices(imgn, posn);
print_convert('tutorial130c.png', '-density 100')
