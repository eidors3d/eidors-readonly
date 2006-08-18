% Explore Stimulation Patterns
% $Id: tutorial030a.m,v 1.1 2006-08-18 17:36:19 aadler Exp $

% 3D Model
imdl_3d= mk_common_model('n3r2',16);
fmdl= imdl_3d.fwd_model;

% Show
fmdl.stimulation=mk_stim_patterns(16,2, '{op}','{op}', ...
             {'meas_current','no_redundant'} );
subplot(121)
display_meas(fmdl,'ya')

fmdl.stimulation=mk_stim_patterns(16,2, '{ad}','{ad}', ...
             {'no_meas_current','do_redundant'} );
subplot(122)
display_meas(fmdl,'ya')

print -r75 -dpng tutorial030a.png;
