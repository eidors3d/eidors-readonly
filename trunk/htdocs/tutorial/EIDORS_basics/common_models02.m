% Distmesh models $Id$
% Distmesh models with fixed electode nodes
imdl=mk_common_model('c2d0d',14);
subplot(221); show_fem(imdl.fwd_model,[0,1]);

imdl=mk_common_model('d2d3d',9);
subplot(222); show_fem(imdl.fwd_model,[0,1]);

print -dpng -r125 common_models02a.png
