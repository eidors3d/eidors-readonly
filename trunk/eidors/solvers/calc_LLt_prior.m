function LLt_prior = calc_LLt_prior( data0, inv_model )
% LLt = calc_LLt_prior( data0, inv_model )
% CALC_LLt_PRIOR: calculate image regularization prior
%   L*L' (which is an estimate of the inverse of the covariance in time)
%
%   Typically, the image prior is matrix n_frames x n_frames of the
%   normalized a priori cross-correlation between measurement frames
%   where frames are a set of measurements taken at nearly the same time
% 
% calc_LLt_prior can be called as
%    LLt_prior= calc_LLt_prior( data0, ... )
%
% and will call the function inv_model.LLt_prior
% parameters to LLt_prior should be passed in the field
% inv_model.LLt_prior_function_name.parameters
%
% If inv_model.LLt_prior is a matrix, calc_LLt_prior will return that matrix,
%
% LLt_prior    the calculated LLt regularization prior
% inv_model    is an inv_model structure
%
% If a function to calculate LLt_prior is not provided,
% LLt = L_prior * L_prior';

% (C) 2017 Alistair Boyle. License: GPL version 2 or version 3
% $Id$

if isfield(inv_model,'LLt_prior')
   if isnumeric(inv_model.LLt_prior)
      LLt_prior = inv_model.LLt_prior;
   else
      try inv_model.LLt_prior = str2func(inv_model.LLt_prior); end
      LLt_prior= eidors_cache( inv_model.LLt_prior, inv_model );
   end
elseif isfield(inv_model,'L_prior')
   LLt_prior = eidors_cache(@calc_from_L_prior, inv_model, 'calc_LLt_prior');
else
   error('calc_LLt_prior: neither L_prior nor LLt_prior provided');
end

function LLt_prior = calc_from_L_prior(inv_model)

   % The user has provided an R prior. We can use this to
   % calculate LLt= L*L';
   if isnumeric(inv_model.L_prior)
      L = inv_model.L_prior;
   else
      try inv_model.L_prior = str2func(inv_model.L_prior); end
      L= eidors_cache( inv_model.L_prior, inv_model );
   end

   LLt_prior = L*L';
