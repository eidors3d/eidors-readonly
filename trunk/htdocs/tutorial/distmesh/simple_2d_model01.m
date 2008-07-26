% Simple distmesh model $Id$

% Create small model with little mesh refinement
imdm=mk_common_model('a2d1c',16);
subplot(121)
show_fem(imdm.fwd_model);
axis on; axis equal

% close in view
subplot(122)
show_fem(imdm.fwd_model);
axis equal; axis([-.1 .1 .85 1.05]);

print -dpng -r125 simple_2d_model01a.png
