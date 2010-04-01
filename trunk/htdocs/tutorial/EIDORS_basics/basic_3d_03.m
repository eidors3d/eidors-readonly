% Basic 3d model $Id$

% Add a circular object at 0.2, 0.5
% Calculate element membership in object

img2.calc_colours.cb_shrink_move = [0.3,0.6, 0];
subplot(221); show_fem(img2,1);

print -dpng -r125 basic_3d_02a.png
