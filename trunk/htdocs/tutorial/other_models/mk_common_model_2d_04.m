% 2D models $Id$

clf; subplot(221)

imdl= mk_common_model('a2d0c',8);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_04a.png

imdl= mk_common_model('b2d1c',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_04b.png

imdl= mk_common_model('c2d2c',32);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_04c.png

imdl= mk_common_model('e2d4c',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_04d.png

