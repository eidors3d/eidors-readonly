% 2D models $Id: mk_common_model_2d_01.m,v 1.1 2008-05-19 17:55:38 aadler Exp $

clf; subplot(221)

imdl= mk_common_model('a2c0',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_01a.png

imdl= mk_common_model('b2c0',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_01b.png

imdl= mk_common_model('c2c0',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_01c.png

imdl= mk_common_model('d2c0',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_01d.png
