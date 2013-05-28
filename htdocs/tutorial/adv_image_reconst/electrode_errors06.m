% GREIT - error $Id$

imdl = eidors_obj('inv_model','','fwd_model',fmdl);
imdl.meas_icov_rm_elecs.elec_list = 13;
imdl.meas_icov_rm_elecs.exponent  = -1;
imdl.meas_icov_rm_elecs.SNR       = 100;
opt.noise_covar = meas_icov_rm_elecs(imdl);
imdl=mk_GREIT_model(img, 0.25, [], opt);

imgr= inv_solve(imdl, zc_resp(:,1), zc_resp(:,22));
imgr.calc_colours.ref_level=0; show_fem(imgr,[0,1]);

axis equal; axis([-1,1.05,-.8,.8]);
print_convert electrode_errors06a.png
