% 2D solver $Id$

imdl.hyperparameter.value= 0.05;

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,(1-.3)*scl];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc0= inv_solve(imdl, vh, vi);
% Show image of reconstruction in upper planes
show_slices(imgc0);
print_convert centre_slice04b.png '-density 75';

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,(1+.3)*scl];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc1= inv_solve(imdl, vh, vi);

% Show image of reconstruction in lower planes
show_slices(imgc1);
print_convert centre_slice04c.png '-density 75';
