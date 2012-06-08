% 2D solver $Id$

% Create a classic 2D inverse model
imdl= mk_common_model('c2c2',16);
% Rotate electrodes to match
imdl.fwd_model.electrode([5:-1:1,16:-1:6]) = imdl.fwd_model.electrode;
imdl.RtR_prior = @prior_gaussian_HPF;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value= 0.01;

imgc= inv_solve(imdl, vh, vi);

% Show on the mesh
subplot(121); show_fem(imgs); axis equal
subplot(122); show_fem(imgc); axis equal
print_convert square_mesh06a.png

% Show on a grid 
subplot(121); show_slices(imgs); axis equal
subplot(122); show_slices(imgc); axis equal
print_convert square_mesh06b.png
