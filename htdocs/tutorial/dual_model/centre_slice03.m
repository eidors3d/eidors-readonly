% 2D solver $Id$


% Set coarse as reconstruction model
imdl.rec_model= c_mdl;
c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,scl];
c_mdl.mk_coarse_fine_mapping.f2c_project = (1/scl)*speye(3);
c_mdl.mk_coarse_fine_mapping.z_depth = inf;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @prior_gaussian_HPF;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value= 0.01;

imgc= inv_solve(imdl, vh, vi);

clf; show_slices(imgc);
print_convert centre_slice04a.png '-density 75';
