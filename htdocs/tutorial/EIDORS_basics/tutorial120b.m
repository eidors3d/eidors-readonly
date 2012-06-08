% Compare 2D algorithms
% $Id$

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;

% This is not an inverse crime; inv_mdl != fwd_mdl
imb=  mk_common_model('b2c',16);
inv2d.fwd_model= imb.fwd_model;
inv2d.fwd_model.np_fwd_solve.perm_sym= '{y}';

% Guass-Newton solvers
inv2d.solve=       @np_inv_solve;

% Tikhonov prior
inv2d.hyperparameter.value = 1e-3;
inv2d.RtR_prior=   @tikhonov_image_prior;
imgr(1)= inv_solve( inv2d, vh, vi);
imgn(1)= inv_solve( inv2d, vh, vi_n);

% NOSER prior
inv2d.hyperparameter.value = 1e-2;
inv2d.RtR_prior=   @noser_image_prior;
imgr(2)= inv_solve( inv2d, vh, vi);
imgn(2)= inv_solve( inv2d, vh, vi_n);

% Laplace image prior
inv2d.hyperparameter.value = 1e-3;
inv2d.RtR_prior=   @laplace_image_prior;
imgr(3)= inv_solve( inv2d, vh, vi);
imgn(3)= inv_solve( inv2d, vh, vi_n);

% Automatic hyperparameter selection
inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'value');
inv2d.hyperparameter.func = @choose_noise_figure;
inv2d.hyperparameter.noise_figure= 0.5;
inv2d.hyperparameter.tgt_elems= 1:4;
inv2d.RtR_prior=   @prior_gaussian_HPF;
inv2d.solve=       @inv_solve_diff_GN_one_step;
imgr(4)= inv_solve( inv2d, vh, vi);
imgn(4)= inv_solve( inv2d, vh, vi_n);
inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'func');

% Total variation using PDIPM
inv2d.hyperparameter.value = 1e-5;
inv2d.solve=       @TV_diffusivity_solve;
inv2d.R_prior=     @prior_TV;
inv2d.parameters.max_iterations= 10;
inv2d.parameters.term_tolerance= 1e-3;

imgr(5)= inv_solve( inv2d, vh, vi);
imgn(5)= inv_solve( inv2d, vh, vi_n);

% Output image
imgn(1).calc_colours.npoints= 128;
imgr(1).calc_colours.npoints= 128;

show_slices(imgr, [inf,inf,0,1,1]);
print_convert tutorial120b.png;

show_slices(imgn, [inf,inf,0,1,1]);
print_convert tutorial120c.png;

