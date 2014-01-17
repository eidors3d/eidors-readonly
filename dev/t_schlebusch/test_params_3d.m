%% create forward model:
fmdl= ng_mk_cyl_models([0.3,0.15,0.01],[16,[0.1,0.2]],[0.01 0 0.002]); 
% since the prior pushes to a ball of zero radius at (0,0,0), we need that
% to be a plausible location
fmdl.nodes(:,3) = fmdl.nodes(:,3) - 0.15;

% Create stimulation pattern:
stim =  mk_stim_patterns(16,2,[0,1],[0,1], {'no_meas_current'}, 1);
fmdl.stimulation = stim;

imdl = select_imdl(fmdl,{'Basic GN abs'});

img = mk_image(fmdl,1,'conductivity');

vh  = fwd_solve(img);

circ= @(r,x,y,z,v) 1 - v * elem_select(fmdl,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',x,y,z,r),'x','y','z'));

bladder = @(r,v) 1 - v * elem_select(fmdl,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',-r-0.01,0,r+0.01,r),'x','y','z'));

bladder_r = 0.03;
b_x = 0.02;
b_y = 0.03;
b_z = 0.03;
% simplified physiological bladder anatomy:
%img.conductivity.elem_data = circ(bladder_r, (0.15-bladder_r-0.01), 0, (bladder_r+0.01), 0.5); % r, x, y, z, cond


% in-plane test model:
% img.conductivity.elem_data = circ(bladder_r, b_x, b_y, b_z, 0.5); % r, x, y, z, cond
img.conductivity.elem_data = bladder(bladder_r,0.5);

vi1  = fwd_solve(img);
vi = add_noise(1e2,vi1);
disp(norm(calc_difference_data(vi1,vi,fmdl)));
show_fem(img)

%% Reconstruction
bkgnd.conductivity.params = [0.01]'; % r x y z
bkgnd.physics_data_mapper = 'conductivity';
bkgnd.conductivity.func   = @apply_bladder_parametrization;
bkgnd.physics_param_mapper = {'conductivity.params'};
imdl.jacobian_bkgnd = bkgnd;
imdl.fwd_model.jacobian = @jacobian_perturb;
imdl.RtR_prior = @prior_tikhonov;
imdl.hyperparameter.value = 1e-9;
imdl.fwd_model.jacobian_perturb.delta  = 0.001;
imdl.parameters.max_iterations = 5;
imdl.parameters.show_iterations  = 1;
imdl.inv_solve.calc_solution_error = 0;
imdl.inv_solve_abs_GN.do_starting_estimate = 0;

rimg = inv_solve(imdl, vi);
figure
show_fem(rimg);
eidors_colourbar(rimg);
disp(rimg.conductivity.params);

