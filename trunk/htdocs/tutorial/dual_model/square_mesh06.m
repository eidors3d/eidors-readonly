% 2D solver $Id: square_mesh06.m,v 1.1 2008-03-28 17:55:23 aadler Exp $

% Create a classic 2D inverse model
imdl= mk_common_model('c2c2',16);
% Rotate electrodes to match
imdl.fwd_model.electrode([5:-1:1,16:-1:6]) = imdl.fwd_model.electrode;
imdl.RtR_prior = @gaussian_HPF_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.01;

imgc= inv_solve(imdl, vh, vi);

% Show on the mesh
subplot(121); show_fem(imgs);
subplot(122); show_fem(imgc);
print -r150 -dpng square_mesh06a.png

% Show on a grid 
subplot(121); show_slices(imgs);
subplot(122); show_slices(imgc);
print -r150 -dpng square_mesh06b.png
