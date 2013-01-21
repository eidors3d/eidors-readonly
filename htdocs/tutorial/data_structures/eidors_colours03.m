% Show EIDORS colours $Id$
subplot(221); img1= img;
img1.calc_colours.cb_shrink_move = [0.5,0.8,.02];

subplot(221)
show_fem(img1,1);
axis equal; axis off; 
print_convert eidors_colours03a.png '-density 75'

subplot(221)
img1.calc_colours.clim= 1;
show_fem(img1,1);
axis equal; axis off;
print_convert eidors_colours03b.png '-density 75'

subplot(221)
img1.calc_colours.clim= 0.3;
show_fem(img1,1);
axis equal; axis off; 
print_convert eidors_colours03c.png '-density 75'
