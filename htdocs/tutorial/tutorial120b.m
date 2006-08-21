% Compare 2D algorithms
% $Id: tutorial120b.m,v 1.2 2006-08-21 04:30:45 aadler Exp $

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;

% This is not an inverse crime; inv_mdl != fwd_mdl
imb=  mk_common_model('b2c',16);
inv2d.fwd_model= imb.fwd_model;
inv2d.fwd_model.misc.perm_sym= '{y}';

% Guass-Newton solver
inv2d.hyperparameter.value = 1e-3;
inv2d.RtR_prior=   @laplace_image_prior;
inv2d.solve=       @np_inv_solve;
imgr(1)= inv_solve( inv2d, vi, vh);
imgn(1)= inv_solve( inv2d, vi_n, vh);

% Automatic hyperparameter selection
inv2d.hyperparameter.func = @aa_calc_noise_figure;
inv2d.hyperparameter.noise_figure= 2;
inv2d.hyperparameter.tgt_elems= 1:4;
inv2d.RtR_prior=   @aa_calc_image_prior;
inv2d.solve=       @aa_inv_solve;
imgr(2)= inv_solve( inv2d, vi, vh);
imgn(2)= inv_solve( inv2d, vi_n, vh);

% Total variation using PDIPM
inv2d.hyperparameter.value = 1e-2;
inv2d.solve=       @ab_tv_diff_solve;
inv2d.R_prior=     @ab_calc_tv_prior;
inv2d.parameters.max_iterations= 20;
inv2d.parameters.term_tolerance= 1e-3;

% TVimg will add the background value
tvimg= inv_solve( inv2d, vi, vh);
tvimg.elem_data = tvimg.elem_data - bkgnd;
imgr(3)= tvimg;

tvimg= inv_solve( inv2d, vi_n, vh);
tvimg.elem_data = tvimg.elem_data - bkgnd;
imgn(3)= tvimg;

% Output image
subplot(211)
show_slices(imgr, [inf,inf,0,1,1]);
subplot(212)
show_slices(imgn, [inf,inf,0,1,1]);

print -r100 -dpng tutorial120b.png;
