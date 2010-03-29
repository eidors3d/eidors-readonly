% Distmesh models $Id$
imdl=mk_common_model('c2d0c',14);
subplot(221); show_fem(imdl.fwd_model,[0,1]);

imdl=mk_common_model('d2d3c',9);
subplot(222); show_fem(imdl.fwd_model,[0,1]);

print -dpng -r125 common_models01a.png
