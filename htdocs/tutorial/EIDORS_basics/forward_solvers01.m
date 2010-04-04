% Forward solvers $Id$

% 2D Model
imdl= mk_common_model('d2d1c',19);

% Create an homogeneous image
img_1 = mk_image(imdl);
subplot(221); show_fem(img_1);

% Add a circular object at 0.2, 0.5
% Calculate element membership in object
img_2 = img_1;
select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
img_2.elem_data = 1 + elem_select(img_2.fwd_model, select_fcn);

img_2.calc_colours.cb_shrink_move = [0.3,0.6,+0.03];
subplot(222); show_fem(img_2,1);

print_convert forward_solvers01a.png
