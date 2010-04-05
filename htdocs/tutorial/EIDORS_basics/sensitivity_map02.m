% Sensitivity map $Id$

imdl= mk_common_model('a2c0',16); % Basic model - replace
imdl.fwd_model= ng_mk_cyl_models([2,1,.07],[16,1.0],[0.1]);
imdl.fwd_model.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
J= calc_jacobian(calc_jacobian_bkgnd(imdl));
img = mk_image(imdl, J(5,:)');

img.calc_colours.clim= 3e-6;

img.calc_colours.npoints= 256;
img.calc_slices.filter = conv2(ones(3),ones(3));
show_3d_slices(img,[0.7,1.0],[0.5],[]);
view(-70,22);
print_convert('sensitivity_map02a.png','-density 60');
view(-70,62);
print_convert('sensitivity_map02b.png','-density 60');

np = interp_mesh(img.fwd_model);
img.elem_data(np(:,3)>1.05)= 0;
show_fem(img);

% Show better
view(-70,52);
crop_model(gca, inline('-x/5 + z>1.05','x','y','z'))
print_convert('sensitivity_map02c.png','-density 60');
