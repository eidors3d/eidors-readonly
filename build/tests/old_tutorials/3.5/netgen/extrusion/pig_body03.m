img_v = img;
% Stimulate between elecs 16 and 5 to get more interesting pattern
img_v.fwd_model.stimulation(1).stim_pattern = sparse([16;5],1,[1,-1],16,1);
img_v.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img_v);

img_v = rmfield(img, 'elem_data');
img_v.node_data = vh.volt(:,1);
img_v.calc_colours.npoints = 128;
img_v.calc_colours.clim = 1.2;

subplot(221);
show_slices(img_v,[inf,inf,0.8]); axis off; axis equal
print_convert pig_body03a.jpg
show_slices(img_v,[inf,inf,1.0]); axis off; axis equal
print_convert pig_body03b.jpg
show_slices(img_v,[inf,inf,1.2]); axis off; axis equal
print_convert pig_body03c.jpg
