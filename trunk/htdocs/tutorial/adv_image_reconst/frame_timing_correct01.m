fmdl = mk_library_model('neonate_16el_lungs');
[fmdl.stimulation fmdl.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}',{'rotate_meas'});
fmdl = mdl_normalize(fmdl, 1);  % Use normalized difference imaging
opt.noise_figure = 0.2; opt.target_size = 0.1;
opt.square_pixels = 1; opt.distr = 3;

inv_mdl = mk_GREIT_model(fmdl, 0.25,5, opt); 
