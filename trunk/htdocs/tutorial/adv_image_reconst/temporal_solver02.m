% Image reconstruction of moving objects $Id$

time_steps=  3;
time_weight= .8;

base_model = mk_common_model( 'c2c2', 16 ); % 576 element
base_model.fwd_model = mdl_normalize(base_model.fwd_model, 0);
% HP value gives good solutions (calculated via Noise Figure)
hp= sqrt(0.101046); % hp is squared
base_model.hyperparameter.value= hp;

% GN Solver
imdl_GN = base_model;
imdl_GN.RtR_prior= @prior_noser;
imdl_GN.prior_noser.exponent= .5;
imdl_GN.solve= @inv_solve_diff_GN_one_step;
imdl_GN.hyperparameter.value= hp;
imdl_GN.fwd_model = mdl_normalize(imdl_GN.fwd_model, 0);

% Temporal Solver
imdl_TS = base_model;
imdl_TS.RtR_prior= @prior_time_smooth;
imdl_TS.prior_time_smooth.space_prior= @prior_noser;
imdl_TS.prior_noser.exponent= .5;
imdl_TS.prior_time_smooth.time_weight= time_weight;
imdl_TS.prior_time_smooth.time_steps=  time_steps;
imdl_TS.solve= @inv_solve_time_prior;
imdl_TS.inv_solve_time_prior.time_steps=   time_steps;

% Kalman Solver
imdl_KS = base_model;
imdl_KS.RtR_prior= @prior_noser;
imdl_KS.prior_noser.exponent= .5;
imdl_KS.solve= @inv_solve_diff_kalman;
