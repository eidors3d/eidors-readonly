% Simulate obj $Id$

% GREIT v1
i_gr = mk_common_gridmdl('GREITc1');

show_fem( inv_solve( i_gr, vh, vt) ); axis equal;
print_convert 'GREIT_examples02a.png' '-density 50'; 

% current GREIT recommendation
opt.noise_figure = 0.5;
opt.distr = 0; % best for cylinders
i_grc = mk_GREIT_model(fmdl,.2,[],opt);
show_fem( inv_solve( i_grc, vh, vt) ); axis equal;
print_convert 'GREIT_examples02b.png' '-density 50'; 

% Sheffield Backprojection
i_bp = mk_common_gridmdl('backproj');

show_fem( inv_solve( i_bp, vh, vt) ); axis equal;
print_convert 'GREIT_examples02c.png' '-density 50';

% 2D Gauss Newton Inverse
i_gn = mk_common_model('d2c2',16);
i_gn.hyperparameter.value = 0.01;
i_gn.fwd_model = mdl_normalize(i_gn.fwd_model, 1);
% i_gn.RtR_prior = @prior_gaussian_HPF;

show_fem( inv_solve( i_gn, vh, vt) ); axis equal;
print_convert 'GREIT_examples02d.png' '-density 50';

% Test the Noise Figure of the GN inverse => 0.5
% i_gn.hyperparameter.tgt_data.meas_t1 = vh;
% i_gn.hyperparameter.tgt_data.meas_t2 = vt;
% calc_noise_figure(i_gn)

