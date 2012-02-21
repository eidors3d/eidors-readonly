function meas_icov = calc_reciproc_error(inv_model, data )
% CALC_RECIPROC_ERROR: CALCULATE RECIPROCITY ERROR MATRIX
% meas_icov = calc_reciproc_error(inv_model, data )
%
% Calculate meas_icov from data
% 
% specify tau as 
%    inv_model.calc_reciproc_error.tau
% as tau increases, effect is less. There is almost
% no effect when tau > 1e-4
%
% also accepts a fwd_model parameter
%
% Reference: Real-time management of faulty electrodes in
%  electrical impedance tomography AE Hartinger, R Guardo,
%  A Adler, H Gagnon. IEEE T BME 2008.

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

try;    tau = inv_model.calc_reciproc_error.tau;
catch;  tau = 2.5e-6;                            end

switch inv_model.type
  case 'inv_model'; fmdl = inv_model.fwd_model;
  case 'fwd_model'; fmdl = inv_model;
  otherwise;        error('calc_reciproc_error: require model input');
end

data = calc_difference_data( 0, data, fmdl);

idx = reciprocity_idx( fmdl );
if any(isnan(idx));
   error('not all meas have reciprocity for this stim_pattern');
end

data = data/max(abs(data));
recerr=  data - data(idx,:);
recerr= mean( abs(recerr), 2);
s2  = exp(-recerr.^2/tau); % eqn 8 from paper

nmeas = length(idx);
meas_icov= spdiags(s2(:),0,nmeas,nmeas);

function oldcode
% Use only with no_rotate
Nel = length(inv_model.fwd_model.electrode);
nframes= size(data,2);
nmeas = size(data,1);
if nmeas ~= Nel^2
   data2= zeros(Nel^2, nframes);
   mselect = inv_model.fwd_model.meas_select;
   data2(mselect,:)= data;
   data2= reshape( data2, Nel,Nel,nframes);
else 
   data2= reshape( data, Nel,Nel,nframes);
end

ndata2= abs(data2)/max(abs(data2(:)));
ndata2r= permute(ndata2, [2,1,3]); % Transpose
recerr= mean( ndata2 - ndata2r, 3);
e2    = recerr.^2; % eqn 8 from paper

s2  = exp(-e2/tau);

if isfield(inv_model.fwd_model,'meas_select')
   mselect = inv_model.fwd_model.meas_select;
   s2= s2(mselect);
   nmeas = length(find(mselect));
end

meas_icov= spdiags(s2(:),0,nmeas,nmeas);
