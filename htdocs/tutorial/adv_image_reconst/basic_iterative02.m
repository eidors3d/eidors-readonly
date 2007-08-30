% Reconstruct images $Id: basic_iterative02.m,v 1.3 2007-08-30 03:58:27 aadler Exp $

% Set reconstruction parameters
imdl_3d.solve= @np_inv_solve;
imdl_3d.RtR_prior= @np_calc_image_prior;
imdl_3d.np_calc_image_prior.parameters= [3 1];
imdl_3d.hyperparameter.value= .01;

imdl_3d.fwd_model.solve=      @np_fwd_solve;
imdl_3d.fwd_model.jacobian=   @np_calc_jacobian;
imdl_3d.fwd_model.system_mat= @np_calc_system_mat;

% Number of iterations and tolerance (defaults)
imdl_3d.parameters.max_iterations = 1;
imdl_3d.parameters.term_tolerance = 1e-3;

%Add 30dB SNR noise to data
noise_level= std(inh_data.meas - homg_data.meas)/10^(30/20);
noise_level=0;
inh_data.meas = inh_data.meas + noise_level* ...
                randn(size(inh_data.meas));

% Reconstruct Images: 1 Iteration
subplot(131)
imdl_3d.parameters.max_iterations = 1;
rec_img= inv_solve(imdl_3d, homg_data, inh_data);
show_slices(rec_img,slice_posn);


% Reconstruct Images: 2 Iterations
subplot(132)
imdl_3d.parameters.max_iterations = 2;
rec_img= inv_solve(imdl_3d, homg_data, inh_data);
show_slices(rec_img,slice_posn);

% Reconstruct Images: 6 Iterations
subplot(133)
imdl_3d.parameters.max_iterations = 6;
rec_img= inv_solve(imdl_3d, homg_data, inh_data);
show_slices(rec_img,slice_posn);


print -r75 -dpng basic_iterative02a.png
