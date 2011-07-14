% Simple distmesh model $Id$

% Create small model with little mesh refinement
imdm=mk_common_model('d2d4c',16);
subplot(121)
show_fem(imdm.fwd_model);
axis on; axis equal; axis([-1.1 1.1 -1.1 1.1]);

% close in view
subplot(122)
show_fem(imdm.fwd_model);
axis equal; axis([-.1 .1 .85 1.05]);

print_convert simple_2d_model02a.png '-density 125'
