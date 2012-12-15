% Reject electrodes + Image $Id$

fmdl= mk_library_model('adult_male_16el_lungs');
% fmdl.electrode = fmdl.electrode([9:16,1:8]);

[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
fmdl.stimulation = stim;
fmdl = mdl_normalize(fmdl, 1);
fmdl.meas_select = msel;
img = mk_image(fmdl, 1);
img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 0.3; % lungs
clear opt; opt.imgsz = [32 32]; opt.square_pixels=1;
opt.noise_figure = .5;
imdl=mk_GREIT_model(img, 0.25, [], opt);
