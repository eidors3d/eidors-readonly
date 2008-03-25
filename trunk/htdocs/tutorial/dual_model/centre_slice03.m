% 2D solver $Id: centre_slice03.m,v 1.1 2008-03-25 19:21:49 aadler Exp $

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,0];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc0= inv_solve(imdl, vh, vi);
imgc0.fwd_model= c_mdl;

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,.2];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc1= inv_solve(imdl, vh, vi);
imgc1.fwd_model= c_mdl;

subplot(132)
show_slices(imgc0);
subplot(133)
show_slices(imgc1);
