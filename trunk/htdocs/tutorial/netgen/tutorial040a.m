% Use Netgen
% $Id: tutorial040a.m,v 1.1 2007-06-15 18:05:51 aadler Exp $

[tank_mdl,centres] = create_tank_mesh_ng( ...
  15,      ... % tank_radius,
  30,      ... % tank_height,
  'R',     ... % CorR,  - circular or rectangular
  4,       ... % log2_electrodes_per_plane,
  2,       ... % no_of_planes,
  10,      ... % first_plane_starts,
  10,      ... % height_between_centres,
  2,       ... % electrode_width, (radius for circ elecs)
  3,       ... % electrode_height,
  'my_mdl',... % fnstem - filename without extension
  0 );         % elec_mesh_density - increase mesh near electrodes

subplot(121)
show_fem(tank_mdl)
view(-28,10)

[tank_mdl2,centres] = create_tank_mesh_ng( ...
  15, 30, 'C', 4, 2, 10, 10, 2, 3, 'my_mdl', 10 );
subplot(122)
show_fem(tank_mdl2)
view(-28,10)

print -r100 -dpng tutorial040a.png;
