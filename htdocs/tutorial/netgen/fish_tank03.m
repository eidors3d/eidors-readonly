img.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img);

img_v = rmfield(img, 'elem_data');
img_v.node_data = vv.volt(:,1);
img_v.calc_colours.npoints = 128;

show_slices(img_v,[inf,inf,fz/2]);
print_convert fish_tank03a.png '-density 60'

% slice in y between fish and target
show_slices(img_v,[inf, mean([fy,0.5]),inf]);
print_convert fish_tank03b.png '-density 60'
