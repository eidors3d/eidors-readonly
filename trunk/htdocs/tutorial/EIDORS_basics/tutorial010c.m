% Reconstruct images
% $Id$

subplot(131)
show_fem(sim_img);

%Add 20dB SNR noise to data
noise_level= std(inh_data.meas - homg_data.meas)/10^(20/20);
inh_data.meas = inh_data.meas + noise_level* ...
                randn(size(inh_data.meas));

%reconstruct 
rec_img= inv_solve(imdl_3d, homg_data, inh_data);

% Show reconstruction as a 3D mesh
subplot(132)
show_fem(rec_img)

subplot(133)
rec_img.calc_colours.npoints = 128;
show_slices(rec_img,[inf,inf,2.0,1,1; ...
                     inf,inf,1.0,1,2]);

print_convert('tutorial010c.png', '-density 100', 0.5);
