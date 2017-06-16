% Explore Stimulation Patterns
% $Id$

% 3D Model
imdl_3d= mk_common_model('n3r2',[16. 2]);
fmdl= imdl_3d.fwd_model;

% Show opposite pattern
fmdl.stimulation=mk_stim_patterns(16,2, '{op}','{op}', ...
             {'meas_current','no_redundant'} );
subplot(121)
show_stim_meas_pattern(fmdl,'ya')

% Show adjacent pattern
fmdl.stimulation=mk_stim_patterns(16,2, '{ad}','{ad}', ...
             {'no_meas_current','do_redundant'} );
subplot(122)
show_stim_meas_pattern(fmdl,'ya')

print_convert tutorial030a.png '-density 75';
