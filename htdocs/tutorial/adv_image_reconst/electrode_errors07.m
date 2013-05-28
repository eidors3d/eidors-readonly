% Basic GREIT reconst $Id$

% use -ve tau to get covariance rather than inv covar
imdl.calc_reciproc_error.tau = -3e-4;
opt.noise_covar = calc_reciproc_error( imdl, zc_resp );

imdl=mk_GREIT_model(img, 0.25, [], opt);

imgr= inv_solve(imdl, zc_resp(:,1), zc_resp(:,22));
imgr.calc_colours.ref_level=0; show_fem(imgr,[0,1]);

axis equal; axis([-1,1.05,-.8,.8]);

print_convert electrode_errors07a.png
