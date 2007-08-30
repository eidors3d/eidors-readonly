% Show EIDORS colours $Id: eidors_vars02.m,v 1.1 2007-08-30 03:30:42 aadler Exp $

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
