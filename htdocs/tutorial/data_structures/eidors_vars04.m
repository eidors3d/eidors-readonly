% Show EIDORS colours $Id: eidors_vars04.m,v 1.1 2007-08-30 03:30:44 aadler Exp $

idx= [51:52, 56:-1:53, 57:60, 64:-1:61, 49,50]; % clockwise index elements 
img.elem_data(idx)= linspace(-2,2,16);
subplot(131)
calc_colours('clim',1.3); % Limit colours to 1,3
show_slices(img);

subplot(132)
calc_colours('clim',[]); % default
show_slices(img);

subplot(133)
img.calc_colours.clim= .3;
show_slices(img);
print -r100 -dpng eidors_vars04a.png
