% Forward solvers $Id$


% 2D Model
imdl_2d= mk_common_model('d2d1c',19);

% Create an homogeneous image
img_2d = calc_jacobian_bkgnd(imdl_2d);
subplot(221); show_fem(img_2d);

% Add a circular object at 0.2, 0.5
% Calculate element membership in object
xctr = 0.2; yctr = 0.5; rad = 0.1;
pts = interp_mesh( img_2d.fwd_model, 4);
dist = (pts(:,1,:)-xctr).^2 + (pts(:,2,:)-yctr).^2;
member = mean( dist < rad^2 ,3);
img_2d.elem_data = 1 + member;

img.calc_colours.cb_shrink_move = [0.5,0.8,0];
subplot(222); show_fem(img_2d,1);

print -dpng -r125 forward_solvers01a.png
