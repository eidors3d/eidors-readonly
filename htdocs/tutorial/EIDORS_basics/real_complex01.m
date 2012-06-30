imdl = mk_common_model('h2d1c',16);

DeltaC = 0+0.1i;
img = mk_image(imdl,1);
vh = fwd_solve(img);
target= mk_c2f_circ_mapping(img.fwd_model, [0.5;0.0;0.1]);
img.elem_data = 1+ DeltaC*target;
vi = fwd_solve(img);
vi = add_noise(5,vi,vh);

subplot(221);
show_fem(img);
img.calc_colours.component = 'imag';
subplot(222);
show_fem(img);


imdl = mk_common_model('b2d1c',16);
imgr = inv_solve(imdl, vh, vi);
imgr.calc_colours.clim = 0.02;
subplot(223);
show_fem(imgr);
imgr.calc_colours.component = 'imag';
subplot(224);
show_fem(imgr);
