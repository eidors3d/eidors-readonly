% Show EIDORS colours $Id: eidors_vars03.m,v 1.1 2007-08-30 03:30:42 aadler Exp $

clf;
% Set square figure and make figure fill the axis
axes('position',[0 0 1 1]);
pp= get(gcf,'paperposition');
set(gcf,'paperposition',[pp(1:3),pp(3)]);

calc_colours('greylev',.001); % black background level
show_slices(img);
print -r20 -dpng eidors_vars03a.png

calc_colours('greylev',.2); % grey background level
show_slices(img);
print -r20 -dpng eidors_vars03b.png

calc_colours('greylev',-.2); %light grey background level
show_slices(img);
print -r20 -dpng eidors_vars03c.png

calc_colours('greylev',-.001); %white background level (default)
show_slices(img);
print -r20 -dpng eidors_vars03d.png


calc_colours('backgnd',[0.2,0.1,0.15]);
show_slices(img);
calc_colours('backgnd',[0.5,0.5,0.15]); %default value
print -r20 -dpng eidors_vars03e.png
set(gcf,'paperposition',pp(1:4));
