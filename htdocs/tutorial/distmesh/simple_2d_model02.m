% Simple distmesh model $Id: simple_2d_model02.m,v 1.1 2008-03-29 01:23:31 aadler Exp $

% Create small model with little mesh refinement
imdm=mk_common_model('d2d4c',16);
subplot(121)
show_fem(imdm.fwd_model);
axis on; axis equal

% close in view
subplot(122)
show_fem(imdm.fwd_model);
axis equal; axis([-.1 .1 .85 1.05]);

print -dpng -r125 simple_2d_model02a.png
