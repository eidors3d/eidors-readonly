%% Train GREIT
opt.imgsz = [64 64]; % 64-by-64 image (yes, we can do that now)
opt.distr = 3; % non-random, uniform
opt.Nsim = 500; % 500 hundred targets to train on, seems enough
opt.target_size = 0.01; %small targets
opt.target_offset = 0;
opt.noise_figure = 0.5; % this is key!
imdl=mk_GREIT_model(img, 0.25, [], opt);

