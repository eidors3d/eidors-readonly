function [r_est, cond_est] = test_params_3d(r_real, cond_real, imdl_greit)
%r_real = 0.03;
%cond_real = 0.7;
%% create forward model:
fmdl= ng_mk_cyl_models([0.3,0.15,0.005],[32,[0.1,0.2]],[0.01 0 0.002]); 
% since the prior pushes to a ball of zero radius at (0,0,0), we need that
% to be a plausible location
fmdl.nodes(:,3) = fmdl.nodes(:,3) - 0.1;
fmdl.nodes(:,1) = fmdl.nodes(:,1) - 0.075;

% Create stimulation pattern:
stim =  mk_stim_patterns(32,2,[0,1],[0,1], {'no_meas_current'}, 1);
fmdl.stimulation = stim;

imdl = select_imdl(fmdl,{'Basic GN abs'});

img = mk_image(fmdl,0.4,'conductivity'); % 0.4 S/m background conductivita (roughly tissue @ 50 kHz)

vh  = fwd_solve(img);

bladder = @(r,v) 1 + v * elem_select(fmdl,inline(sprintf('(x-%f).^2 + (y-%f).^2 + (z-%f).^2 <= %f^2',0,0,r,r),'x','y','z'));

img.conductivity.elem_data = bladder(r_real, cond_real);

%
vi1  = fwd_solve(img);
vi = add_noise(1e2,vi1);
disp(norm(calc_difference_data(vi1,vi,fmdl)));
%show_fem(img)

%% Reconstruction
bkgnd.conductivity.params = [0.05, 1]'; % r x y z
bkgnd.data_mapper = 'conductivity';
bkgnd.conductivity.func   = @apply_bladder_parametrization;
bkgnd.params_mapper = {'conductivity.params'};
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

r_est = rimg.conductivity.params(1);
cond_est = rimg.conductivity.params(2);
end