% Show EIDORS colours $Id$
subplot(221);
img1= img;
img1.calc_colours.cb_shrink_move = [0.5,0.8,.02];

show_fem(img1,1);
axis equal; axis off; 
print_convert eidors_colours02a.png '-density 75'

img1.calc_colours.greylev= -0.35;
show_fem(img1,1);
axis equal; axis off; 
print_convert eidors_colours02b.png '-density 75'

img1.calc_colours.greylev= 0.35;
show_fem(img1,1);
axis equal; axis off; 
print_convert eidors_colours02c.png '-density 75'
