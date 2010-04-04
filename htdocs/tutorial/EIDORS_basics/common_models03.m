% Distmesh models $Id$
imdl=mk_common_model('c2C1',16);
subplot(221); show_fem(imdl.fwd_model,[0,1]);

imdl=mk_common_model('d2t3',16);
subplot(222); show_fem(imdl.fwd_model,[0,1]);

print_convert common_models03a.png
