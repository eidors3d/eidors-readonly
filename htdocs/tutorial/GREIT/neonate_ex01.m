% Inverse model
%fmdl = mk_library_model('neonate_16el_lungs');
elec_pos = [16,1,.5]; elec_shape=[0.05,0,0,0,60]; maxsz=0.08; nfft=27;
elec_pos = [16,1,.5]; elec_shape=[0.05,0,0,0,60]; maxsz=0.08; nfft=27;
%fmdl = mk_library_model({'neonate','boundary','left_lung','right_lung'},elec_pos, elec_shape, maxsz,nfft);
 fmdl = mk_library_model({'neonate','boundary'},elec_pos, elec_shape, maxsz,nfft);

[fmdl.stimulation fmdl.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl = mdl_normalize(fmdl,1);

img = mk_image(fmdl,1); img.elem_data(vertcat(fmdl.mat_idx{2:3})) = 0.3;
img.calc_colours.ref_level=1;

show_fem(img);
return
opt.square_pixels = 1;
opt.noise_figure = 0.5; opt.imgsz = [64 64];
imdl = mk_GREIT_model(img, 0.20, [], opt);

