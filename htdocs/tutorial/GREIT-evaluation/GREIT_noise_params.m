function params = GREIT_noise_params(imgs, alg, vh, vi )
% params = GREIT_noise_params(imgs, homg_voltage, sig_voltage)
%  params(1,:) = Noise Figure

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% There are better ways here
noise = 0.01*std(vh)*randn(208,1000);
snr_y = mean(abs(vi-vh*ones(1,size(vi,2))),1) / mean(std(noise,[],1),2); 

im_n= feval(alg, vh, vh*ones(1,size(noise,2)) + noise);
snr_x = mean(mean(abs(imgs),1),2) / mean(abs(im_n(:)));
params= [snr_y(:)./snr_x(:)]';

