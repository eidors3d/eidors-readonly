% Show EIDORS colours $Id: eidors_colours05.m,v 1.1 2008-03-14 15:26:02 aadler Exp $
subplot(221); img1= img;

show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours05a.png

img1.calc_colours.window_range= 0.90;
show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours05b.png

img1.calc_colours.window_range= 0.20;
show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours05c.png
