% 2D solver $Id: centre_slice03.m,v 1.2 2008-03-27 19:47:24 aadler Exp $

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,-.3];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc0= inv_solve(imdl, vh, vi);
% Show image of reconstruction in upper planes
subplot(132)
show_slices(imgc0);

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,.3];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc1= inv_solve(imdl, vh, vi);

% Show image of reconstruction in lower planes
subplot(133)
show_slices(imgc1);

print -r150 -dpng centre_slice03a.png;
