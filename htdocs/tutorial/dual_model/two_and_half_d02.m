% Solve 2D and 3D model $Id: two_and_half_d02.m,v 1.4 2008-03-19 19:02:45 aadler Exp $

% Original target
subplot(141)
show_fem(inhomg_img); view(-62,28)

% Create inverse Model: Classic
imdl= eidors_obj('inv_model','NPs GN inverse');
imdl.solve=       @np_inv_solve;
imdl.hyperparameter.value = 1e-3;
imdl.R_prior= @np_calc_image_prior;
imdl.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
imdl.jacobian_bkgnd.value= 1;
imdl.reconst_type= 'difference';
imdl.fwd_model= f_mdl; % fine model

% Classic (inverse crime) solver
img1= inv_solve(imdl, vh, vi);
subplot(142)
show_fem(img1); view(-62,28)
