imdl = eidors_obj('inv_solve','test');

imdl.hyperparameter.value = .3;
if 0
imdl.RtR_prior = @prior_laplace;
imdl.solve = @inv_solve_diff_GN_one_step;
end
imdl.R_prior = @prior_TV;
imdl.solve = @inv_solve_diff_pdipm;
imdl.inv_solve_diff_pdipm.norm_image = 1;
imdl.inv_solve_diff_pdipm.norm_data  = 1;

imdl.jacobian_bkgnd.value = 1;
imdl.reconst_type = 'difference';
imdl.fwd_model = gallery_3D_img.fwd_model;
imdl.rec_model = gallery_3D_img.rec_model;

vh = ref_data.meas;
vi = real_data.meas;
k = vh\vi;

imgr= inv_solve(imdl, vh, vi/k);
imgr.fwd_model = imdl.fwd_model;
