% Sensitivity map $Id$

imdl= mk_common_model('f2d1c',16); % 'j2d1c'
J= calc_jacobian(calc_jacobian_bkgnd(imdl));
img = mk_image(imdl, J(5,:)');

img.calc_colours.npoints= 256;
img.calc_slices.filter = conv2(ones(3),ones(3));
img.calc_colours.clim = .0004;
show_slices(img);
print_convert('sensitivity_map01a.png','-density 60');

img.calc_colours.cb_shrink_move = [0.5,0.8,0];
show_fem(img,[1,1]);
print_convert('sensitivity_map01b.png','-density 60');
