% Show EIDORS colours $Id$
clf; subplot(221); img1= img;
img1.calc_colours.cb_shrink_move = [0.5,0.8,-.10];

clf; subplot(221)
show_fem(img1,1);
axis equal; axis off; axis tight;
print_convert eidors_colours03a.png '-density 75'

clf; subplot(221)
img1.calc_colours.clim= 1;
show_fem(img1,1);
axis equal; axis off; axis tight;
print_convert eidors_colours03b.png '-density 75'

clf; subplot(221)
img1.calc_colours.clim= 0.3;
show_fem(img1,1);
axis equal; axis off; axis tight;
print_convert eidors_colours03c.png '-density 75'
