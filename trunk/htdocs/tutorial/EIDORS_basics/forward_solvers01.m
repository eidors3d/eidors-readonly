% Forward solvers $Id$

% 2D Model
imdl= mk_common_model('d2d1c',19);

% Create an homogeneous image
img_1 = calc_jacobian_bkgnd(imdl);
subplot(221); show_fem(img_1);

% Add a circular object at 0.2, 0.5
% Calculate element membership in object
img_2 = img_1;
xctr = 0.2; yctr = 0.5; rad = 0.1;
pts = interp_mesh( img_2.fwd_model, 4);
dist = (pts(:,1,:)-xctr).^2 + (pts(:,2,:)-yctr).^2;
member = mean( dist < rad^2 ,3);
img_2.elem_data = 1 + member;

img_2.calc_colours.cb_shrink_move = [0.3,0.6,+0.03];
subplot(222); show_fem(img_2,1);

print -dpng -r125 forward_solvers01a.png
