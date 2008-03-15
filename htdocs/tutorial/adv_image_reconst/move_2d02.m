% Generate eidors planar finite element model
mdl2dim = mk_common_model('b2c');
% Set eidors_obj hyperparameter member.
mdl2dim.hyperparameter.value= 0.01;

clim= .088;

% Solve inverse problem for mdl2dim eidors_obj model.
img2dim = inv_solve(mdl2dim, vh, vi);
img2dim.calc_colours.clim= clim;
img2dim.calc_colours.backgnd= [.9,.9,.9];

% Plot results for each algorithm
subplot(1,2,1);
show_fem_move(img2dim);
calc_colours(img2dim,[],1); % do colourbar

% Set eidors_obj hyperparameter member.
mdlM = mdl2dim;
% Place traditional jacobian in temporary member.
mdlM.fwd_model.conductivity_jacobian = mdlM.fwd_model.jacobian;
% Redefine jacobian member for movement & conductivity.
mdlM.fwd_model.jacobian = 'aa_e_move_jacobian';
mdlM.RtR_prior =     'aa_e_move_image_prior';
mdlM.aa_e_move_image_prior.parameters = sqrt(1e2/1); 

% Solve inverse problem for mdlM eidors_obj model.
imgM = inv_solve(mdlM, vh, vi);
imgM.calc_colours.clim= clim;
imgM.calc_colours.backgnd= [.9,.9,.9];

% Plot results for each algorithm
subplot(1,2,2);
show_fem_move(imgM);

calc_colours(imgM,[],1); % do colourbar

set(gcf,'paperposition',[.25 2.5 8 4]);
print -r125 -dpng move_2d02.png
