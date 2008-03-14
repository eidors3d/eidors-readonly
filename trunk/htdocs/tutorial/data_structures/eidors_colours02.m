% Show EIDORS colours $Id: eidors_colours02.m,v 1.1 2008-03-14 15:26:02 aadler Exp $
subplot(221); img1= img;

show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours02a.png

img1.calc_colours.greylev= -0.35;
show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours02b.png

img1.calc_colours.greylev= 0.35;
show_fem(img1,1);
axis equal; axis off; print -r75 -dpng eidors_colours02c.png
