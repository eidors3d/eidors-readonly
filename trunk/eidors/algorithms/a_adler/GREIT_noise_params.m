function params = GREIT_noise_params(imdl, vh, vi )
% params = GREIT_noise_params(imdl, homg_voltage, sig_voltage)
%  params(1,:) = Noise Figure

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% NOTE THAT WE ASSUME A LINEAR ALGORITHM FOR THIS MEASURE

imgs = inv_solve(imdl,vh,vi);
imgs = imgs.elem_data;
% There are better ways here
noise = 0.01*std(vh)*randn(208,1000);
vhn= vh*ones(1,size(noise,2));
%signal_y = vi -  (vh*ones(1,size(vi,2)));
 signal_y = vi ./ (vh*ones(1,size(vi,2))) - 1;
%noise_y  = mean(std(noise      ),2); 
 noise_y  = mean(std(noise./vhn ),2); 
snr_y = mean(abs(signal_y),1) / noise_y;

im_n= inv_solve(imdl, vh, vhn + noise);
im_n = im_n.elem_data;
%signal_x = mean(mean(abs(imgs),1),2);
 signal_x = mean(mean(   (imgs),1),2);
%noise_x  = mean(abs(im_n(:)));
 noise_x  = mean(std(std(im_n)));
snr_x = signal_x / noise_x;
params= [snr_y(:)./snr_x(:)]';
