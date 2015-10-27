img = mk_image(fmdl,1); img.elem_data(vertcat(fmdl.mat_idx{2:3})) = 0.6;
opt.square_pixels = 1; opt.imgsz = [64 64];
%opt.noise_figure = 0.5;
%imdl = mk_GREIT_model(img, 0.20, [], opt);
imdl = mk_GREIT_model(img, 0.20, 10, opt);
