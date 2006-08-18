% Explore Stimulation Patterns
% $Id: tutorial030b.m,v 1.1 2006-08-18 18:14:16 aadler Exp $

% We have a 16 electrode EIT machine with adjacent drive
adjdrv= mk_stim_patterns(16,1, '{ad}','{ad}', ...
             {'no_meas_current','do_redundant'} );

% Arrange 16 electrodes in a zigzag
zigzag_mdl= fmdl;
zigzag_mdl.stimulation= adjdrv;
zigzag_pat= [ 1:2:15;
             18:2:32]; 
zigzag_mdl.electrode= fmdl.electrode( zigzag_pat(:) );

subplot(121);
display_meas(zigzag_mdl,'ya')

% Arrange 16 electrodes as planar
planar_mdl= fmdl;
planar_mdl.stimulation= adjdrv;
planar_pat= [ 1:2:15;
             17:2:32]'; 
planar_mdl.electrode= fmdl.electrode( planar_pat(:) );

subplot(122);
display_meas(planar_mdl,'ya')

print -r75 -dpng tutorial030b.png;
