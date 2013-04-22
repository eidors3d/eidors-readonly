% Create fwd models
% $Id$

subplot(121);

% 2D Model
imdl_2d= mk_common_model('b2c',16);
show_fem(imdl_2d.fwd_model);

axis square
subplot(122);

% 3D Model
imdl_3d= mk_common_model('n3r2',[16,2]);
show_fem(imdl_3d.fwd_model);

axis square; view(-35,14);
print_convert('tutorial010a.png','-density 100')
