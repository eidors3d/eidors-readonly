function params = GREIT_noise_params(imgs, alg, vh, vi )
% params = GREIT_noise_params(imgs, homg_voltage, sig_voltage)
%  params(1,:) = Noise Figure

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% There are better ways here
noise = 0.01*std(vh)*randn(208,1000);
vhn= vh*ones(1,size(noise,2));
%signal_y = vi -  (vh*ones(1,size(vi,2)));
 signal_y = vi ./ (vh*ones(1,size(vi,2))) - 1;
%noise_y  = mean(std(noise      ),2); 
 noise_y  = mean(std(noise./vhn ),2); 
snr_y = mean(abs(signal_y),1) / noise_y;

im_n= feval(alg, vh, vhn + noise);
%signal_x = mean(mean(abs(imgs),1),2);
 signal_x = mean(mean(   (imgs),1),2);
%noise_x  = mean(abs(im_n(:)));
 noise_x  = mean(std(std(im_n)));
snr_x = signal_x / noise_x;
params= [snr_y(:)./snr_x(:)]';
