% Distmesh models $Id$
imdl=mk_common_model('c2d0c',14);
show_fem(imdl.fwd_model,[0,1]);

print_convert('common_models01a.png','-density 60');

imdl=mk_common_model('d2d3c',9);
show_fem(imdl.fwd_model,[0,1]);

print_convert('common_models01b.png','-density 60');
