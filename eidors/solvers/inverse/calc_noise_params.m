function params = calc_noise_params(imdl, vh, vi )
% params = GREIT_noise_params(imdl, homg_voltage, sig_voltage)
%  params(1,:) = Noise Figure = SNR(image) / SNR(data)
%
%  see also: eval_GREIT_fig_merit or using test_performance

% (C) 2008 Andy Adler. Licensed under GPL v2 or v3
% $Id$

% NOTE THAT WE ASSUME A LINEAR ALGORITHM FOR THIS MEASURE

if ischar(imdl) && strcmp(imdl,'UNIT_TEST'), do_unit_test, return, end

if 0 % old code with random noise
     % Keep this to validate while we test it
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

try 
    VOL = get_elem_volume(imdl.rec_model);
catch
    VOL = get_elem_volume(imdl.fwd_model);
end
VOL = spdiags(VOL,0, length(VOL), length(VOL));

signal_x = VOL*signal_x;
noise_x = VOL*noise_x;

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

function do_unit_test
ll = eidors_msg('log_level',1);
% test1; % cannot deal with c2f
imdl = mk_common_model('a2t2',16); test2(imdl);
imdl = mk_common_model('d2t2',16); test2(imdl);
imdl = mk_common_model('a2c2',16); test2(imdl);
imdl = mk_common_model('d2c2',16); test2(imdl);
ll = eidors_msg('log_level',ll);

function test1
% big model with c2f and supposedly an NF of 0.5
fmdl = mk_library_model('pig_23kg_16el');
[fmdl.stimulation fmdl.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl = mdl_normalize(fmdl, 1);  % Use normalized difference imaging
opt.noise_figure = 0.5; opt.imgsz = [64 64];
imdl = mk_GREIT_model(fmdl, 0.25, [], opt);
% homogeneous measurement
img = mk_image(fmdl,1);
vh = fwd_solve(img);
% inhomogeneous measurement
select_fcn = inline('(x-0).^2+(y-0).^2+(z-0.5).^2<0.1^2','x','y','z');
mfrac = elem_select(fmdl, select_fcn);
img.elem_data = img.elem_data + mfrac*0.1;
vi = fwd_solve(img);

nf1 = calc_noise_params(imdl, vh.meas, vi.meas);

imdl.hyperparameter.tgt_data.meas_t1 = vh.meas;
imdl.hyperparameter.tgt_data.meas_t2 = vi.meas;
try
    % calc_noise_figure doens't support dual models
    nf2 = calc_noise_figure(imdl);
catch
    nf2 = 0;
end
unit_test_cmp('Noise fig implementations',nf1, nf2, 1e-2);

function test2(imdl)
fmdl = imdl.fwd_model;
% homogeneous measurement
img = mk_image(fmdl,1);
vh = fwd_solve(img);
% inhomogeneous measurement
select_fcn = inline('(x-0).^2+(y-0).^2.^2<15^2','x','y','z');
mfrac = elem_select(fmdl, select_fcn);
img.elem_data = img.elem_data + mfrac*0.1;
vi = fwd_solve(img);

nf1 = calc_noise_params(imdl, vh.meas, vi.meas);
eidors_msg(nf1,0);

imdl.hyperparameter.tgt_data.meas_t1 = vh.meas;
imdl.hyperparameter.tgt_data.meas_t2 = vi.meas;
% calc_noise_figure doens't support dual models
nf2 = calc_noise_figure(imdl,[],1000);
unit_test_cmp('Noise fig implementations',nf1, nf2, 1e-2);
