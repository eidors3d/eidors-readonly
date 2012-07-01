% $Id$
opt.imgsz = [32 32];
opt.distr = 3; % non-random, uniform
opt.Nsim = 1000;
opt.target_size = 0.05; % Target size (frac of medium)
opt.noise_figure = 1.0; % Recommended NF=0.5;

fmdl = ng_mk_ellip_models([2, 2,1.4,0.2] ,[n_elecs,layers],[0.1]);
fmdl.stimulation =  stim;
fmdl = mdl_normalize(fmdl, 0);
img = mk_image(fmdl,1);

imdl = mk_GREIT_model(img, 0.25, [], opt);

show_fem(fmdl); view(0,70);

print_convert mk_GREIT_mat_2layer02a.png '-density 60'
