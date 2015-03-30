% Generate eidors planar finite element model
mdl2dim = mk_common_model('b2c');
mdl2dim.hyperparameter.value= 0.1;

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

mdlM.fwd_model.jacobian = @jacobian_movement;
mdlM.RtR_prior =          @prior_movement;
mdlM.prior_movement.parameters = sqrt(1e2/1); 

% Solve inverse problem for mdlM eidors_obj model.
imgM = inv_solve(mdlM, vh, vi);
imgM.calc_colours.clim= clim;
imgM.calc_colours.backgnd= [.9,.9,.9];

% Plot results for each algorithm
subplot(1,2,2);
show_fem_move(imgM);
imgM.calc_colours.cb_shrink_move = [0.5,0.9,.02];
eidors_colourbar(imgM);

print_convert move_2d02.png
