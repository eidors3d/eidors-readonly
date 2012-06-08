% Reconstruct 3D movement $Id$

img3dim= img;

mdl3dim.RtR_prior = 'prior_laplace';
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

% This is also available by this code:
% mdlM = select_imdl(mdl3dim, {'Elec Move GN'});

% Inverse solution of data with movement consideration
move_vs_conduct = 20;  % Movement penalty (symbol mu in paper)

mdlM = mdl3dim;
if 0
   mdlM.fwd_model.conductivity_jacobian = mdlM.fwd_model.jacobian;
   mdlM.fwd_model.jacobian = @jacobian_movement_perturb; % perturbation type jacobian
else
   mdlM.fwd_model.jacobian = @jacobian_movement;
end

mdlM.RtR_prior = @prior_movement;

mdlM.prior_movement.parameters = move_vs_conduct;

% Solve inversglobale problem and show slices
imgM = inv_solve(mdlM, vh, vi);
imgM.calc_colours.backgnd=[.8 .8 .9];
subplot(1,3,3)
show_slices_move( imgM );
print_convert  move_3d02.png '-density 100'
