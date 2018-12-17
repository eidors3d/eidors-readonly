% Model of 32x1 electrode belt
fmdl= ng_mk_ellip_models([4,0.8,1.1,.5],[32,2.0],[0.05]);
% Swisstom BBVet stimulation pattern
skip4 = {32,1,[0,5],[0,5],{'no_meas_current_next1'},1};
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip4{:});
