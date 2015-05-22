%% Forward model
fmdl= ng_mk_cyl_models([3,2,.4],[16,1,2],[.1,0,.025]);
fmdl.stimulation = mk_stim_patterns(16,2,'{ad}','{ad}');