% Basic GREIT reconst $Id$

opt.noise_covar = 1./(meas_icov_rm_elecs( fmdl,13)+.001);
imdl=mk_GREIT_model(img, 0.25, [], opt);
imdl.fwd_model.meas_select = msel;

imgr= inv_solve(imdl, zc_resp(:,1), zc_resp(:,22));
imgr.calc_colours.ref_level=0; show_fem(imgr,[0,1]);

axis equal; axis([-1,1.05,-.8,.8]);
print_convert electrode_errors06a.jpg
