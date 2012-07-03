% Distmesh models $Id$
% Distmesh models with fixed electode nodes
imdl=mk_common_model('c2d0d',14);
show_fem(imdl.fwd_model,[0,1]);

print_convert('common_models02a.png','-density 60');

imdl=mk_common_model('d2d3d',9);
show_fem(imdl.fwd_model,[0,1]);

print_convert('common_models02b.png','-density 60');
