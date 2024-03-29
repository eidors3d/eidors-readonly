[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
img.fwd_model.stimulation = stim;
img.fwd_model.normalize_measurements=1;
opt.imgsz = [32 32];
opt.distr = 3; % non-random, uniform
opt.Nsim = 500; % 500 hundred targets to train on, seems enough
opt.target_size = 0.03; %small targets
opt.target_offset = 0;
opt.noise_figure = .5; % this is key!
imdl=mk_GREIT_model(img, 0.25, [], opt);
imdl.fwd_model.meas_select = msel;
