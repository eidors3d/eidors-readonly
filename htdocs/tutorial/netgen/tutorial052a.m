% Radial move ball in a single electrode plane
% $Id: tutorial052a.m,v 1.1 2007-06-15 18:05:51 aadler Exp $

electrodes_per_plane= 16;
number_of_planes= 1;
movement_pattern= 'radial_move';
finelevel= '-veryfine';
refine_electrodes= 20;
tank_radius= 15;
tank_height= 4;
electrode_width = 0.25;
electrode_height= 0.25;
rect_or_circ_electrode= 'C';

move_the_ball;
