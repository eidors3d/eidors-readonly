% Basic GREIT reconst $Id$

img= inv_solve(imdl, zc_resp(:,1), zc_resp(:,22));
img.calc_colours.ref_level=0;
show_slices(img);

opt.noise_covar = 1./(meas_icov_rm_elecs( fmdl,13)+.001);
imdl=mk_GREIT_model(img, 0.25, [], opt);
imdl.fwd_model.meas_select = msel;

print_convert electrode_errors05a.jpg
