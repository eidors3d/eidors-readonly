% Inverse model
fmdl = mk_library_model('neonate_16el');
[fmdl.stimulation fmdl.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl = mdl_normalize(fmdl,1);
img = mk_image(fmdl,1);
opt.noise_figure = 0.5; opt.imgsz = [64 64];
imdl = mk_GREIT_model(img, 0.25, [], opt);
% imdl = mk_common_gridmdl('GREITc1');

% Data: eidors3d.sf.net/data_contrib/if-neonate-spontaneous/if-neonate-spontaneous.zip
vv= eidors_readdata('P04P-1016.get');
