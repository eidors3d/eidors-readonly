% Reconstruct images
% $Id: tutorial010c.m,v 1.1 2007-06-15 18:17:51 aadler Exp $

subplot(131)
show_fem(sim_img);

%Add 20dB SNR noise to data
noise_level= std(inh_data.meas - homg_data.meas)/10^(20/20);
inh_data.meas = inh_data.meas + noise_level* ...
                randn(size(inh_data.meas));

%reconstruct 
rec_img= inv_solve(imdl_3d, inh_data, homg_data);

% Show reconstruction as a 3D mesh
subplot(132)
show_fem(rec_img)

subplot(133)
show_slices(rec_img,[inf,inf,2.0,1,1; ...
                     inf,inf,1.0,1,2]);

pp=get(gcf,'paperposition');
set(gcf,'paperposition',pp.*[1,1,1.4,1]);
print -r75 -dpng tutorial010c.png
set(gcf,'paperposition',pp.*[1,1,1,1]);
