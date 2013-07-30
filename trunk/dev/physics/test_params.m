
fmdl= ng_mk_cyl_models([0,1,0.1],32,[1/6,0,0.03]);
fmdl.stimulation = mk_stim_patterns(32,1,[0 14],[0 14]);

imdl = select_imdl(fmdl,{'Basic GN abs'});

img = mk_image(fmdl,1,'conductivity');

vh  = fwd_solve(img);

circ= @(r,x,y,v) 1 - v * elem_select(fmdl,inline(sprintf('(x-%f).^2 + (y-%f).^2 <= %f^2',x,y,r),'x','y','z'));

img.conductivity.elem_data = circ(0.2, 0.3, 0.7, 0.5);
vi  = fwd_solve(img);
vi = add_noise(20,vi);

bkgnd.conductivity.params = [0.5, 0, 0]'; % r x y
bkgnd.physics_data_mapper = 'conductivity';
bkgnd.conductivity.func   = @apply_circle_parametrization;
bkgnd.physics_param_mapper = {'conductivity.params'};
imdl.jacobian_bkgnd = bkgnd;
imdl.fwd_model.jacobian = @jacobian_perturb;
imdl.RtR_prior = @prior_tikhonov;
imdl.hyperparameter.value = 1e-9;
imdl.fwd_model.jacobian_perturb.delta  = 0.1;
imdl.parameters.max_iterations = 10 ;
imdl.parameters.show_iterations  = 1;

rimg = inv_solve(imdl, vi);
show_fem(rimg);
eidors_colourbar(rimg);
disp(rimg.conductivity.params);

