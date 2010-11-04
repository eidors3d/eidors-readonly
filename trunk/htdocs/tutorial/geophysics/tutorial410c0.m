imdl = eidors_obj('inv_solve','test');

imdl.hyperparameter.value = .1;
if 0
imdl.RtR_prior = @laplace_image_prior;
imdl.solve = @aa_inv_solve;
end
imdl.R_prior = @ab_calc_tv_prior;
imdl.solve = @pdipm_diff;
imdl.pdipm_diff.norm_image = 1;
imdl.pdipm_diff.norm_data  = 1;

imdl.jacobian_bkgnd.value = 1;
imdl.reconst_type = 'difference';
imdl.fwd_model = gallery_3D_img.fwd_model;
imdl.rec_model = gallery_3D_img.rec_model;

vh = ref_data.meas;
vi = real_data.meas;
k = vh\vi;

imgr= inv_solve(imdl, vh, vi/k);
imgr.fwd_model = imdl.fwd_model;
