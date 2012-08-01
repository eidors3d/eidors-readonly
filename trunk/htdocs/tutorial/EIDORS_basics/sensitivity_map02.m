% Sensitivity map $Id$

fwd_model= ng_mk_cyl_models([2,1,.07],[16,1.0],[0.1]);
fwd_model.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
J= calc_jacobian( mk_image(fwd_model,1) );
Sens = J(5,:)'./get_elem_volume(fwd_model);
img = mk_image(fwd_model, Sens');

img.calc_colours.clim= 3e-2;

img.calc_colours.npoints= 256;
% img.calc_slices.filter = conv2(ones(3),ones(3));
img.calc_colours.transparency_thresh = -1;
show_3d_slices(img,[0.7,1.0],[0.5],[]);
view(-70,22);
print_convert('sensitivity_map02a.png','-density 60');
view(-70,62);
print_convert('sensitivity_map02b.png','-density 60');

np = interp_mesh(img.fwd_model);
img.calc_colours.transparency_thresh = .15;
img.elem_data(np(:,3)>1.05)= 0;
show_fem(img);

% Show better
view(-70,52);
crop_model(gca, inline('-x/5 + z>1.05','x','y','z'))
print_convert('sensitivity_map02c.png','-density 60');
