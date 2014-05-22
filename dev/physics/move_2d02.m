% Generate eidors planar finite element model
mdl2dim = mk_common_model('b2c');
mdl2dim.hyperparameter.value= 0.01;
mdl2dim.inv_solve_diff_GN_one_step.calc_step_size = 1;
% mdl2dim.inv_solve_diff_GN_one_step.step_size = 1;

clim= .088;

% Solve inverse problem for mdl2dim eidors_obj model.
img2dim = inv_solve(mdl2dim, vh, vi);
img2dim.calc_colours.clim= clim;
img2dim.calc_colours.backgnd= [.9,.9,.9];

% Plot results for each algorithm
subplot(1,2,1);
show_fem_move(img2dim);
img2dim.calc_colours.cb_shrink_move = [0.5,1.0,.02];
eidors_colourbar(img2dim);

% Set eidors_obj hyperparameter member.
mdlM = mdl2dim;

mdlM.fwd_model.solve    = @fwd_solve_elec_move;
mdlM.fwd_model.jacobian = @jacobian_movement;
mdlM.RtR_prior =          @prior_movement;
mdlM.prior_movement.parameters = sqrt(1e2/1); 
mdlM.jacobian_bkgnd = rmfield(mdlM.jacobian_bkgnd, 'value');
mdlM.jacobian_bkgnd.conductivity.elem_data = ones(length(mdlM.fwd_model.elems),1);
mdlM.jacobian_bkgnd.movement.electrode_data = zeros(32,1);
mdlM.jacobian_bkgnd.data_mapper = 'conductivity';
mdlM.jacobian_bkgnd.params_mapper = {'conductivity.elem_data' 'movement.electrode_data'};

% Solve inverse problem for mdlM eidors_obj model.
imgM = inv_solve(mdlM, vh, vi);
imgM.calc_colours.clim= clim;
imgM.calc_colours.backgnd= [.9,.9,.9];

% Plot results for each algorithm
subplot(1,2,2);
show_fem(imgM);
imgM.calc_colours.cb_shrink_move = [0.5,0.9,.02];
eidors_colourbar(imgM);

print_convert move_2d02.png
