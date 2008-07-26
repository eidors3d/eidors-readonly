% Reconstruct 3D movement $Id$

img3dim= img;

mdl3dim.RtR_prior = 'laplace_image_prior';
mdl3dim.hyperparameter.value = 3e-3;

% Show slices of 3D model with true movement vectors
subplot(1,3,1);
img3dim.elem_data = img3dim.elem_data - 1;
img3dim.calc_colours.backgnd=[.8 .8 .9];
show_slices_move( img3dim, move );

% Inverse solution of data without movement consideration
img3dim= inv_solve(mdl3dim, vh, vi);
img3dim.calc_colours.backgnd=[.8 .8 .9];
subplot(1,3,2)
show_slices_move( img3dim );

% Inverse solution of data with movement consideration
move_vs_conduct = 20;  % Movement penalty (symbol mu in paper)
% Define a eidglobalors_obj Movement model solved by electrode movement 
% algorithms.
mdlM = mdl3dim;
mdlM.fwd_model.conductivity_jacobian = mdlM.fwd_model.jacobian;

% this is a perturbation type jacobian
%mdlM.fwd_model.jacobian = 'aa_e_move_jacobian';

mdlM.fwd_model.jacobian = @calc_move_jacobian; % faster / more accurate
mdlM.RtR_prior = 'aa_e_move_image_prior';
mdlM.aa_e_move_image_prior.parameters = move_vs_conduct;
% Solve inversglobale problem and show slices
imgM = inv_solve(mdlM, vh, vi);
imgM.calc_colours.backgnd=[.8 .8 .9];
subplot(1,3,3)
show_slices_move( imgM );
print -r100 -dpng move_3d02.png
