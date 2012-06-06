function params = calc_noise_params(imdl, vh, vi )
% params = GREIT_noise_params(imdl, homg_voltage, sig_voltage)
%  params(1,:) = Noise Figure = SNR(image) / SNR(data)
%
%  see also: eval_GREIT_fig_merit or using test_performance

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% NOTE THAT WE ASSUME A LINEAR ALGORITHM FOR THIS MEASURE

if 0 % old code with random noise
   Nnoise = 1000;
   noise = 0.01*std(vh)*randn(size(vh,1),Nnoise);
   vhn= vh*ones(1,Nnoise) + noise;
else % use independent noise model on each channel
   noise = 0.01*std(vh)*speye(size(vh,1));
   vhn= vh*ones(1,size(vh,1)) + noise;
end

signal_y = calc_difference_data( vh, vi,  imdl.fwd_model);
noise_y  = calc_difference_data( vh, vhn, imdl.fwd_model);

signal_x = inv_solve(imdl, vh, vi);  signal_x = signal_x.elem_data;
noise_x  = inv_solve(imdl, vh, vhn); noise_x  = noise_x.elem_data;

signal_x = mean(abs(signal_x),1);
noise_x  = mean(std(noise_x));
snr_x = signal_x / noise_x;

signal_y = mean(abs(signal_y),1);
noise_y  = mean(std(noise_y)); 
snr_y = signal_y / noise_y;

params= [snr_y(:)./snr_x(:)]';

% NOTES on the calculations: AA - Feb 20, 2012
% SNR = mean(abs(x)); VAR = 
% ym= E[y]                
% Sy= E[(y-ym)*(y-ym)'] = E[y*y'] - ym*ym'
% ny = sqrt(trace(Sy))
% xm= E[x]  = E[R*y] = R*E[y] = R*ym
% Sx= E[(x-xm)*(x-xm)'] = E[x*x'] - xm*xm'
%   = E[R*ym*ym'*R'] = R*E[ym*ym']*R' = R*Sy*R'
% nx = sqrt(trace(Sx))
% 
% signal = mean(abs(x));
% 
% In this case, these are exactly the same:
%    noise_x  = mean(std(noise_x));
%    noise_x  = sqrt(mean(noise_x.^2,2));

