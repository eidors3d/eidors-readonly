% Solve 2D and 3D model $Id: two_and_half_d02.m,v 1.2 2008-03-19 00:04:28 aadler Exp $

% Original target
subplot(131)
show_fem(demo_img); view(-62,28)

% Create inverse Model: Classic
imdl.type=       'inv_model';
imdl.name=       'Nick Polydorides GN inverse';
imdl.solve=       @np_inv_solve;
imdl.hyperparameter.value = 3e-1;
imdl.R_prior= @np_calc_image_prior;
imdl.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
imdl.jacobian_bkgnd.value= 1;
imdl.reconst_type= 'difference';
imdl.fwd_model= demo_img.fwd_model;

% Classic (inverse crime) solver
img1= inv_solve(imdl, vh, vi);
subplot(132)
show_fem(img1); view(-62,28)

print -r125 -dpng two_and_half_d02a.png
