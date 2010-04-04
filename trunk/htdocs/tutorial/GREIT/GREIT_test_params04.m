% GREIT params eval $Id$

% Reconstruct GREIT Images
i_bp = mk_common_gridmdl('backproj');
imgr = inv_solve(i_bp, vh, vi);

imgr.calc_colours.npoints = 32;
params = eval_GREIT_fig_merit(imgr, xyzr);

plot(r, params(1,:)); axis([0,0.9,0,2.0]);  ylabel('AR');
print_convert('GREIT_test_params04_ar.png',  '-density 60',0.3);

plot(r, params(2,:)); axis([0,0.9,-0.15,0.15]);  ylabel('PE');
print_convert('GREIT_test_params04_pe.png',  '-density 60',0.3);

plot(r, params(3,:)); axis([0,0.9,0,0.4]);  ylabel('RES');
print_convert('GREIT_test_params04_res.png', '-density 60',0.3);

plot(r, params(4,:)); axis([0,0.9,0,0.3]);  ylabel('SD');
print_convert('GREIT_test_params04_sd.png',  '-density 60',0.3);

plot(r, params(5,:)); axis([0,0.9,0,0.6]);  ylabel('RNG');
print_convert('GREIT_test_params04_rng.png', '-density 60',0.3);
