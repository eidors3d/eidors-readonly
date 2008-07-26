% Show EIDORS colours $Id$

subplot(131)
calc_colours('npoints',32);
show_slices(img);

subplot(132)
calc_colours('npoints',128);
show_slices(img);

subplot(133)
calc_colours('npoints',64);
show_slices(img); %default value
print -r100 -dpng eidors_vars02a.png
