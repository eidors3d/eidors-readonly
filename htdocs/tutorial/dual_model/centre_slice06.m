% 2D solver $Id$

% Set coarse as reconstruction model
imdl.rec_model= c_mdl;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @prior_gaussian_HPF;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value= 0.2;

imgc= inv_solve(imdl, vh, vi);
clf; show_slices(imgc);
print_convert centre_slice06a.png '-density 75';

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,-.3];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc0= inv_solve(imdl, vh, vi);

imgc0= inv_solve(imdl, vh, vi);
show_slices(imgc0);
print_convert centre_slice06b.png '-density 75';

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,.3];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;

imgc1= inv_solve(imdl, vh, vi);
show_slices(imgc1);
print_convert centre_slice06c.png '-density 75';
