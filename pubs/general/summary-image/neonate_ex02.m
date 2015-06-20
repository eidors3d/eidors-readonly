img = mk_image(fmdl,1); img.elem_data(vertcat(fmdl.mat_idx{2:3})) = 0.3;
opt.square_pixels = 1; opt.imgsz = [128 128];
 opt.noise_figure = 0.5;
%imdl = mk_GREIT_model(img, 0.20, [], opt);
 imdl = mk_GREIT_model(img, 0.20, 15, opt);
