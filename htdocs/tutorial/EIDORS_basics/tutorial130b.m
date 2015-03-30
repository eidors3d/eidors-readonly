% Compare 3D algorithms
% $Id$
clear imgr imgn

% Create Inverse Model
inv3d= eidors_obj('inv_model', 'EIT inverse');
inv3d.reconst_type= 'difference';
inv3d.jacobian_bkgnd.value = 1;
inv3d.fwd_model= imb.fwd_model;
inv3d.hyperparameter.value = 0.03; 


% Gauss-Newton Solver
inv3d.solve=       @inv_solve_diff_GN_one_step;

% Tikhonov prior
inv3d.R_prior=     @prior_tikhonov;
imgr(1)= inv_solve( inv3d, vh, vi);
imgn(1)= inv_solve( inv3d, vh, vi_n);

% Laplace prior
inv3d.R_prior=     @prior_laplace;
% inv3d.np_calc_image_prior.parameters= [3 1]; %  deg=1, w=1
imgr(2)= inv_solve( inv3d, vh, vi);
imgn(2)= inv_solve( inv3d, vh, vi_n);

% Andrea Borsic's PDIPM TV solver
inv3d.prior_TV.alpha2 = 1e-5;
inv3d.parameters.max_iterations= 20;
inv3d.parameters.term_tolerance= 1e-3;
inv3d.R_prior=     @prior_TV;
inv3d.solve=       @inv_solve_TV_pdipm;

imgr(3)= inv_solve( inv3d, vh, vi);
imgn(3)= inv_solve( inv3d, vh, vi_n);

% Output image
posn= [inf,inf,2.5,1,1;inf,inf,1.5,1,2;inf,inf,0.5,1,3];
clf;
imgr(1).calc_colours.npoints= 128;
show_slices(imgr, posn);
print_convert('tutorial130b.png', '-density 100')

imgn(1).calc_colours.npoints= 128;
show_slices(imgn, posn);
print_convert('tutorial130c.png', '-density 100')
