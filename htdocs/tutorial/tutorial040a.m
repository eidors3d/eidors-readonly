% Use Netgen
% $Id: tutorial040a.m,v 1.1 2006-08-18 19:33:25 aadler Exp $

[tank_mdl,centres] = create_tank_mesh_ng( ...
  15,     ... % tank_radius,
  30,     ... % tank_height,
  'C',    ... % CorR,  - circular or rectangular
  4,      ... % log2_electrodes_per_plane,
  2,      ... % no_of_planes,
  10,     ... % first_plane_starts,
  10,     ... % height_between_centres,
  2,      ... % electrode_width, (radius for circ elecs)
  3,      ... % electrode_height,
  'my_mdl');  % fnstem - filename without extension


print -r75 -dpng tutorial040a.png;
