% Lung images
% $Id: tutorial310a.m,v 1.1 2006-11-13 18:03:37 aadler Exp $

% 2D Model
imdl= mk_common_model('c2t2',16);

% Make correct stimulation pattern
[st, els]= mk_stim_patterns(...
   16, ... % electrodes / ring
    1, ... % 1 ring of electrodes
   '{ad}','{ad}', ... % adj stim and measurement
   { 'no_meas_current', ... %  don't mesure on current electrodes
     'no_rotate_meas',  ... %  don't rotate meas with stimulation
     'do_redundant', ...    %  do redundant measurements
   }, 10 );  % stimulation current (mA)
imdl.fwd_model.stimulation= st;
imdl.fwd_model.meas_select= els;

% most EIT systems image best with normalized difference
imdl.fwd_model.normalized_measurements= 1;


subplot(221);
show_fem(imdl.fwd_model);

axis equal
print -r100 -dpng tutorial310a.png;
