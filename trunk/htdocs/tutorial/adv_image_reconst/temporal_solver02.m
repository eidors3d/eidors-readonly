% Image reconstruction of moving objects $Id: temporal_solver02.m,v 1.3 2007-08-30 03:58:27 aadler Exp $

time_steps=  3;
time_weight= .8;

base_model = mk_common_model( 'c2c2', 16 ); % 576 element
base_model.fwd_model.normalize_measurements= 0;
% HP value gives good solutions (calculated via Noise Figure)
hp= 0.101046;
base_model.hyperparameter.value= hp;

% GN Solver
imdl_GN = base_model;
imdl_GN.RtR_prior= @noser_image_prior;
imdl_GN.noser_image_prior.exponent= .5;
imdl_GN.solve= @np_inv_solve;
imdl_GN.hyperparameter.value= hp;
imdl_GN.fwd_model.normalize_measurements= 0;

% Temporal Solver
imdl_TS = base_model;
imdl_TS.RtR_prior= @time_smooth_prior;
imdl_TS.time_smooth_prior.space_prior= @noser_image_prior;
imdl_TS.noser_image_prior.exponent= .5;
imdl_TS.time_smooth_prior.time_weight= time_weight;
imdl_TS.time_smooth_prior.time_steps=  time_steps;
imdl_TS.solve= @time_prior_solve;
imdl_TS.time_prior_solve.time_steps=   time_steps;

% Kalman Solver
imdl_KS = base_model;
imdl_KS.RtR_prior= @noser_image_prior;
imdl_KS.noser_image_prior.exponent= .5;
imdl_KS.solve= @inv_kalman_diff;
