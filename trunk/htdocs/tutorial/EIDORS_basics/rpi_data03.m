% RPI tank model $Id$

% simple inverse model -> replace fields to match this model
imdl.fwd_model = fmdl; 

img = inv_solve(imdl , vh, vi);
img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img.calc_colours.ref_level = 0;
show_fem(img,[1,1]); axis off; axis image

print_convert('rpi_data03a.png','-density 60');
