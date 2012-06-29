% Voltage distribution $Id$

nodes = img.fwd_model.nodes; mn = min(nodes); mx= max(nodes);
img.fwd_model.mdl_slice_mapper.x_pts = linspace(mn(1),mx(1),50);
img.fwd_model.mdl_slice_mapper.y_pts = linspace(mn(2),mx(2),10);
img.calc_colours.clim = 5;

show_fem(img); axis image;
hold on;
show_current( img); 
hold off;

print_convert conformal1_03a.png '-density 125'

nodes = img2.fwd_model.nodes; mn = min(nodes); mx= max(nodes);
img2.fwd_model.mdl_slice_mapper.x_pts = linspace(mn(1),mx(1),50);
img2.fwd_model.mdl_slice_mapper.y_pts = linspace(mn(2),mx(2),50);
img2.calc_colours.clim = 5;

show_fem(img2); axis image;
hold on;
show_current( img2); 
hold off;

print_convert conformal1_03a.png '-density 125'
