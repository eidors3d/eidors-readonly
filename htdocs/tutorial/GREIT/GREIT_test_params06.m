% GREIT params eval $Id$

% Reconstruct GREIT Images
i_gr = mk_common_gridmdl('GREITc1');
imgr = inv_solve(i_gr, vh, vi);

imgr.calc_colours.npoints = 32;
params = eval_GREIT_fig_merit(imgr, xyzr);

plot(r, params(1,:)); axis([0,1,0,1.5]);  ylabel('AR');
print_convert('GREIT_test_params06a_ar.png','',0.4);

plot(r, params(2,:)); axis([0,1,-0.2,0.5]);  ylabel('PE');
print_convert('GREIT_test_params06a_pe.png','',0.4);

plot(r, params(3,:)); axis([0,1,0,0.6]);  ylabel('RES');
print_convert('GREIT_test_params06a_res.png','',0.4);

plot(r, params(4,:)); axis([0,1,0,0.6]);  ylabel('SD');
print_convert('GREIT_test_params06a_sd.png','',0.4);

plot(r, params(5,:)); axis([0,1,0,0.6]);  ylabel('RNG');
print_convert('GREIT_test_params06a_rng.png','',0.4);

i_bp = mk_common_gridmdl('backproj');
imgr = inv_solve(i_bp, vh, vi);

imgr.calc_colours.npoints = 32;
params = eval_GREIT_fig_merit(imgr, xyzr);

plot(r, params(1,:)); axis([0,1,0,1.5]);  ylabel('AR');
print_convert('GREIT_test_params06b_ar.png','',0.4);

plot(r, params(2,:)); axis([0,1,-0.2,0.5]);  ylabel('PE');
print_convert('GREIT_test_params06b_pe.png','',0.4);

plot(r, params(3,:)); axis([0,1,0,0.6]);  ylabel('RES');
print_convert('GREIT_test_params06b_res.png','',0.4);

plot(r, params(4,:)); axis([0,1,0,0.6]);  ylabel('SD');
print_convert('GREIT_test_params06b_sd.png','',0.4);

plot(r, params(5,:)); axis([0,1,0,0.6]);  ylabel('RNG');
print_convert('GREIT_test_params06b_rng.png','',0.4);
