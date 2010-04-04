% Simulate obj $Id$

% GREIT v1
i_gr = mk_common_gridmdl('GREITc1');

show_fem( inv_solve( i_gr, vh, vt) );
print_convert 'GREIT_examples02a.png';

% Sheffield Backprojection
i_bp = mk_common_gridmdl('backproj');

show_fem( inv_solve( i_bp, vh, vt) );
print_convert 'GREIT_examples02b.png';

% 2D Gauss Newton Inverse
i_gn = mk_common_model('d2c2',16);
i_gn.hyperparameter.value = 0.1;
i_gn.fwd_model.normalize_measurements = 1;
% i_gn.RtR_prior = @gaussian_HPF_prior;

show_fem( inv_solve( i_gn, vh, vt) );
print_convert 'GREIT_examples02c.png';

% Test the Noise Figure of the GN inverse => 0.5
% i_gn.hyperparameter.tgt_data.meas_t1 = vh;
% i_gn.hyperparameter.tgt_data.meas_t2 = vt;
% calc_noise_figure(i_gn)

