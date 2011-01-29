function vv = add_noise( SNR, v1, v2, options)
% ADD_NOISE: Add a given SNR to EIDORS data
% v1_w_noise = add_noise( SNR, v1, v2, options)
%
% Usage:
%  add_noise( SNR, v1 )           - add noise to v1 where signal = vv
%  add_noise( SNR, v1, v2)        - add noise to v1 where signal = v1 - v2
%  add_noise( SNR, v1, v2,'norm') - add noise to v1 where signal = (v1-v2)/v2
%
% SNR is defined in terms of power SNR =  || signal || / || noise ||

% (C) 2010 Andy Adler. License: GPL v2 or v3. $Id$

if isstr(SNR) && strcmp(SNR,'UNIT_TEST'); do_unit_test; return; end

if nargin>=2; try; v1 = v1.meas; end; end
if nargin>=3; try; v2 = v2.meas; end; end
   
if nargin == 2;
   signal = v1;
elseif nargin==3;
   signal = v1 - v2;
elseif strcmp( options, 'norm' )
   signal = (v1 - v2) ./ v1;
else
   error('add_noise: input arguments not understood');
end
   
noise = randn(size(signal));

% SNR = norm(signal)/norm(noise)
% so  scale norm(noise) by norm(signal)/SNR

noise = noise *norm(signal) / norm(noise) / SNR;

vv = eidors_obj('data','from add_noise');
vv.meas = v1 + noise;

function do_unit_test
    v1 = 2.0*ones(10,1);
    v2 = 2.1*ones(10,1);
    
    v0 = add_noise( 2, v1);
    SNR_test(2, v0.meas - v1, v1);

    v0 = add_noise(.1, v1, v2);
    SNR_test(.1, v0.meas - v1, v1 - v2);

    v0 = add_noise(.3, v1, v2, 'norm' );
    SNR_test(.3, v0.meas - v1, (v1 - v2)./v1);

function SNR_test(SNRspec, noi, sig)
    SNR = norm(sig)/norm(noi);
    if abs(SNR - SNRspec) < .001; disp('ok');
    else;                         disp('fail'); end
