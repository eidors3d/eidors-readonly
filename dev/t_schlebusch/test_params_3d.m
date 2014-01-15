%% create forward model:
fmdl= ng_mk_cyl_models([0.3,0.15,0.01],[16,[0.1,0.2]],[0.01]); 
% Create stimulation pattern:
stim =  mk_stim_patterns(16,2,[0,1],[0,1], {'no_meas_current'}, 1);
fmdl.stimulation = stim;

imdl = select_imdl(fmdl,{'Basic GN abs'});

img = mk_image(fmdl,1,'conductivity');

vh  = fwd_solve(img);

circ= @(r,x,y,z,v) 1 - v * elem_select(fmdl,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',x,y,z,r),'x','y','z'));

bladder_r = 0.02;

% simplified physiological bladder anatomy:
%img.conductivity.elem_data = circ(bladder_r, (0.15-bladder_r-0.01), 0, (bladder_r+0.01), 0.5); % r, x, y, z, cond

% in-plane test model:
img.conductivity.elem_data = circ(bladder_r, 0, 0, 0.15, 0.5); % r, x, y, z, cond

vi  = fwd_solve(img);
vi = add_noise(40,vi);

show_fem(img)

%% Reconstruction
bkgnd.conductivity.params = [0.05, 0, 0, 0.15]'; % r x y z
bkgnd.physics_data_mapper = 'conductivity';
bkgnd.conductivity.func   = @apply_circle_parametrization_3d;
bkgnd.physics_param_mapper = {'conductivity.params'};
imdl.jacobian_bkgnd = bkgnd;
imdl.fwd_model.jacobian = @jacobian_perturb;
imdl.RtR_prior = @prior_tikhonov;
imdl.hyperparameter.value = 1e-9;
imdl.fwd_model.jacobian_perturb.delta  = 0.03;
imdl.parameters.max_iterations = 50 ;
imdl.parameters.show_iterations  = 1;

rimg = inv_solve(imdl, vi);
show_fem(rimg);
eidors_colourbar(rimg);
disp(rimg.conductivity.params);

