% Reconstruct images $Id$

% Set reconstruction inv_solve_gn
imdl_3d.solve= @inv_solve_gn;
imdl_3d.RtR_prior= @prior_laplace;
imdl_3d.hyperparameter.value= .01;
imdl_3d.inv_solve_gn.return_working_variables = 1;
iter_res = @(img) [size(img.inv_solve_gn.r,1)-1, img.inv_solve_gn.r(end,1)];


% Number of iterations and tolerance (defaults)
imdl_3d.inv_solve_gn.max_iterations = 1;

%Add 30dB SNR noise to data
noise_level= std(inh_data.meas - homg_data.meas)/10^(30/20);
noise_level=0;
inh_data.meas = inh_data.meas + noise_level* ...
                randn(size(inh_data.meas));

% Reconstruct Images: 1 Iteration
subplot(131)
imdl_3d.inv_solve_gn.max_iterations = 1;
rec_img= inv_solve(imdl_3d, homg_data, inh_data);
show_slices(rec_img,slice_posn);
title(sprintf('iter=%d resid=%5.3f',iter_res(rec_img)));


% Reconstruct Images: 2 Iterations
subplot(132)
imdl_3d.inv_solve_gn.max_iterations = 2;
rec_img= inv_solve(imdl_3d, homg_data, inh_data);
show_slices(rec_img,slice_posn);
title(sprintf('iter=%d resid=%5.3f',iter_res(rec_img)));

% Reconstruct Images: 5 Iterations -- but stops at 4 (not improving)
subplot(133)
imdl_3d.inv_solve_gn.max_iterations = 5;
rec_img= inv_solve(imdl_3d, homg_data, inh_data);
show_slices(rec_img,slice_posn);
title(sprintf('iter=%d resid=%5.3f',iter_res(rec_img)));


print -r75 -dpng basic_iterative02a.png
