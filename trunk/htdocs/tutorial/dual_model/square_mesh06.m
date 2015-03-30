% 2D solver $Id$

% Create a classic 2D inverse model
imdl= mk_common_model('c2c2',16);
imdl.RtR_prior = @prior_gaussian_HPF;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value= 0.1;

imgc= inv_solve(imdl, vh, vi);

% Show on the mesh
subplot(121); show_fem(imgs); axis equal; axis tight;
subplot(122); show_fem(imgc); axis equal; axis tight;
print_convert square_mesh06a.png

% Show on a grid 
subplot(121); show_slices(imgs)
subplot(122); show_slices(imgc)
print_convert square_mesh06b.png
