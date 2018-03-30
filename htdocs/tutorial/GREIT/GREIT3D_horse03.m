% 2D GREIT model
clear opt;
opt.imgsz = [32 32];
opt.square_pixels = true;
opt.noise_figure = 0.5;
img = mk_image(fmdl,1);
imdl2b= mk_GREIT_model(img, 0.25, [], opt);
