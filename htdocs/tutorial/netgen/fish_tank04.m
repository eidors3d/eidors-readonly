% zoom in on the voltage in a cut plane on the face
img_v = rmfield(img, 'elem_data');
img_v.node_data = vv.volt(:,1);
img_v.calc_colours.npoints = 1000;

show_slices(img_v,[fl/2-0.1,inf,inf]);
axis([700,1000,400,600])
print_convert fish_tank04a.png '-density 60'



