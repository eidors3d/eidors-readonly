% Solve 2D and 3D model $Id: two_and_half_d02.m,v 1.1 2008-03-15 19:13:57 aadler Exp $

% Create inverse Model: Classic
imdl.name= 'Nick Polydorides GN inverse';
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


% Create inverse Model: Classic
imdl.name= 'Coarse/Fine inverse';
imdl.coarse_fine.solve = imdl.solve;
imdl.solve = @coarse_fine_solve;
% imdl.coarse_fine.mapping = speye(n_e3d); % original solver
  imdl.coarse_fine.mapping = c2f;
imdl= eidors_obj('inv_model', imdl);

img2= inv_solve(imdl, vh, vi);
subplot(133)
show_fem(img2); view(-62,28)

% Original
subplot(131)
show_fem(demo_img); view(-62,28)

print -r125 -dpng two_and_half_d02a.png
