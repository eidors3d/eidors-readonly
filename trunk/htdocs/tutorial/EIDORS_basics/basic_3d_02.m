% Basic 3d model $Id$

% Add a circular object at 0.2, 0.5
% Calculate element membership in object
select_fcn = inline('(x-0.2).^2 + (y-0.5).^2 + (z-2).^2 < 0.3^2','x','y','z');
memb_frac = elem_select( img1.fwd_model, select_fcn);
img2 = mk_image(img1, 1 + memb_frac );

img2.calc_colours.cb_shrink_move = [0.3,0.6, 0];
show_fem(img2,1);

print_convert('basic_3d_02a.png','-density 60');
