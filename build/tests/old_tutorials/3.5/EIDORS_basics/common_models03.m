% Distmesh models $Id$
imdl=mk_common_model('c2C1',16);
show_fem(imdl.fwd_model,[0,1]);

print_convert('common_models03a.png','-density 60');

imdl=mk_common_model('d2t3',16);
show_fem(imdl.fwd_model,[0,1]);

print_convert('common_models03b.png','-density 60');
