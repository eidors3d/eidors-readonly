% Element Solvers Target $Id$

% Coarse model
imdl1= mk_common_model('a2c0',16);
imdl1.RtR_prior = @laplace_image_prior;
imdl1.hyperparameter.value = 0.01;

img1 = inv_solve(imdl1, vh, vi);
show_fem(img1);
print_convert nodal_solve02a.png

% Fine model
imdl2= mk_common_model('c2c0',16);
imdl2.RtR_prior = @laplace_image_prior;
imdl2.hyperparameter.value = 0.01;

img2 = inv_solve(imdl2, vh, vi);
show_fem(img2);
print_convert nodal_solve02b.png
