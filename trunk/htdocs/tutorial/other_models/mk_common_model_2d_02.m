% 2D models $Id: mk_common_model_2d_02.m,v 1.1 2008-05-19 17:55:38 aadler Exp $

clf; subplot(221)

imdl= mk_common_model('c2c0',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_02a.png

imdl= mk_common_model('c2C0',16);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_02b.png

imdl= mk_common_model('c2c2',8);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_02c.png

imdl= mk_common_model('c2C2',8);
show_fem(imdl.fwd_model);
print -dpng -r100 mk_common_model_2d_02d.png
