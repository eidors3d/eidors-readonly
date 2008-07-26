% Show EIDORS colours $Id$
subplot(221); img1= img;

show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours04a.png

img1.calc_colours.sat_adj= 0.99;
show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours04b.png

img1.calc_colours.sat_adj= 0.8;
show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours04c.png
