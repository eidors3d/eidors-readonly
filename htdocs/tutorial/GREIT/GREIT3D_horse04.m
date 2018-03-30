vopt.imgsz = [32 32];
vopt.square_pixels = true;
vopt.zvec = linspace(-1,1,10)*1.2+2;
vopt.save_memory = 1;
opt.noise_figure = 1.0;

% GREIT 3D with 2x16 electrode belt
[imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
imdl3= mk_GREIT_model(imdl, 0.20, [], opt);
