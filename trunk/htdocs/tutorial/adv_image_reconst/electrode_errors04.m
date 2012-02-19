% Reject electrodes + Image $Id$

fmdl= mk_library_model('adult_male_16el_lungs');
fmdl.electrode = fmdl.electrode([9:16,1:8]);

[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
fmdl.stimulation = stim;
fmdl.normalize_measurements=1;
img = mk_image(fmdl, 1);
img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 0.3; % lungs

opt.imgsz = [32 32];
opt.square_pixels=1;
opt.noise_figure = .5;
opt.noise_covar = 1./(meas_icov_rm_elecs( fmdl,13)+.001);
imdl=mk_GREIT_model(img, 0.25, [], opt);
imdl.fwd_model.meas_select = msel;

load montreal_data_1995
img= inv_solve(imdl, zc_resp(:,1), zc_resp(:,22));
img.calc_colours.ref_level=0;
show_slices(img);
