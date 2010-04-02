% Element Solvers Target $Id$

% Coarse model
imdl1= mk_common_model('a2c0',16);
imdl1.solve = @nodal_solve;
imdl1.RtR_prior = @laplace_image_prior;
imdl1.hyperparameter.value = 0.01;

img1 = inv_solve(imdl1, vh, vi);
show_fem(img1);
print -dpng -r125 nodal_solve03a.png

% Fine model
imdl2= mk_common_model('c2c0',16);
imdl2.solve = @nodal_solve;
imdl2.RtR_prior = @laplace_image_prior;
imdl2.hyperparameter.value = 0.01;

img2 = inv_solve(imdl2, vh, vi);
show_fem(img2);
print -dpng -r125 nodal_solve03b.png
