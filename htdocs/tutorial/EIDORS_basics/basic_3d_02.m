% Basic 3d model $Id$

% Add a circular object at 0.2, 0.5
% Calculate element membership in object
img2 = img1;
xctr = 0.2; yctr = 0.5; zctr = 2; rad = 0.1;
pts = interp_mesh( img2.fwd_model, 3);
dist = (pts(:,1,:)-xctr).^2 + (pts(:,2,:)-yctr).^2 + (pts(:,3,:)-zctr).^2;
member = mean( dist < rad^2 ,3);
img2.elem_data = 1 + member;

img2.calc_colours.cb_shrink_move = [0.3,0.6, 0];
subplot(221); show_fem(img2,1);

print -dpng -r125 basic_3d_02a.png
