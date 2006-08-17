% tutorial1_create_fwd_model
% $Id: tutorial010a.m,v 1.1 2006-08-17 17:57:33 aadler Exp $

subplot(121);

% 2D Model
imdl_2d= mk_common_model('b2c',16);
show_fem(imdl_2d.fwd_model);

axis square
subplot(122);

% 3D Model
imdl_3d= mk_common_model('n3r2',16);
show_fem(imdl_3d.fwd_model);

axis square; view(-35,14);
print -r75 -dpng tutorial1_create_fwd_model.png
