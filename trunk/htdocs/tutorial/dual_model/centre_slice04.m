% 2D solver $Id: centre_slice04.m,v 1.3 2008-03-28 03:12:59 aadler Exp $

imdl.hyperparameter.value= 0.003;

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,(1-.3)*scl];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc0= inv_solve(imdl, vh, vi);
% Show image of reconstruction in upper planes
subplot(132)
show_slices(imgc0);

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,(1+.3)*scl];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc1= inv_solve(imdl, vh, vi);

% Show image of reconstruction in lower planes
subplot(133)
show_slices(imgc1);

print -r150 -dpng centre_slice04a.png;
