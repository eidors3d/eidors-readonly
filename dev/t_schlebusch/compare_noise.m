function [r_real, cond_real, noise, r_est_parametric, cond_est_parametric, gi] = compare_noise(r_real, cond_real, imdl_greit, noise)
%run('C:\DATEN\eidors_3.7.1\startup.m')
%path(path, 'C:\DATEN\eidors\dev\physics');
cur = cd;
cd C:\DATEN\eidors\eidors
eidors_startup({'physics', 't_schlebusch'});
eidors_cache('cache_size',4e9);
cd(cur)
%r_real = 0.03;
%cond_real = 0.7;
%% create forward model:
fmdl= ng_mk_cyl_models([0.3,0.15,0.01],[32,[0.1,0.2]],[0.01 0 0.002]); 
stim =  mk_stim_patterns(32,2,[0,1],[0,1], {'no_meas_current'}, 1);
% since the prior pushes to a ball of zero radius at (0,0,0), we need that
% to be a plausible location
fmdl.nodes(:,3) = fmdl.nodes(:,3) - 0.1;
fmdl.nodes(:,1) = fmdl.nodes(:,1) - 0.075;

% Create stimulation pattern:
fmdl.stimulation = stim;

imdl = select_imdl(fmdl,{'Basic GN abs'});

img = mk_image(fmdl,1,'conductivity');

vh  = fwd_solve(img);

bladder = @(r,v) 1 + v * elem_select(fmdl,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',0,0,r,r),'x','y','z'));

img.conductivity.elem_data = bladder(r_real, cond_real);

%
vi1  = fwd_solve(img);
vi = add_noise(noise,vi1);
%vi=vi1; %test without noise
disp(norm(calc_difference_data(vi1,vi,fmdl)));
%show_fem(img)

%% Reconstruction
bkgnd.conductivity.params = [0.05, 1]'; % r x y z
bkgnd.physics_data_mapper = 'conductivity';
bkgnd.conductivity.func   = @apply_bladder_parametrization;
bkgnd.physics_param_mapper = {'conductivity.params'};
imdl.jacobian_bkgnd = bkgnd;
imdl.fwd_model.jacobian = @jacobian_perturb;
imdl.RtR_prior = @prior_tikhonov;
%imdl.RtR_prior = @(x) spdiag([1 4]);
imdl.hyperparameter.value = 1e-9;
imdl.fwd_model.jacobian_perturb.delta  = 0.001;
imdl.parameters.max_iterations = 5;
imdl.parameters.show_iterations  = 1;
imdl.inv_solve.calc_solution_error = 0;
imdl.inv_solve_abs_GN.do_starting_estimate = 0;

rimg = inv_solve(imdl, vi);
%figure
%show_fem(rimg);
%eidors_colourbar(rimg);
disp(rimg.conductivity.params);

r_est_parametric = rimg.conductivity.params(1);
cond_est_parametric = rimg.conductivity.params(2);

%% GI Reconstruction
imdl_greit.inv_solve.calc_solution_error = 0;
rimg_gi = inv_solve(imdl_greit, vh, vi); % reconstruct using forward results and calc'd rec. matrix
gi = sum(sum(rimg_gi.elem_data));


end