function RtR_prior = calc_RtR_prior( inv_model, varargin )
% RtR = calc_RtR_prior( inv_model, varargin )
% CALC_RtR_PRIOR: calculate image regularization prior
%   R'*R (which is an estimate of the inverse of the covariance)
%
%   Typically, the image prior is matrix n_elem x n_elem of the
%   normalized a priori crosscorrelation FEM element values
% 
% calc_RtR_prior can be called as
%    RtR_prior= calc_RtR_prior( inv_model, ... )
%
% and will call the function inv_model.RtR_prior
% parameters to RtR_prior should be passed in the field
% inv_model.RtR_prior_function_name.parameters
%
%
% RtR_prior    the calculated RtR regularization prior
% inv_model    is an inv_model structure
%
% If a function to calculate RtR_prior is not provided,
% RtR = R_prior' * R_prior;

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_RtR_prior.m,v 1.5 2007-08-29 09:00:55 aadler Exp $

RtR_prior = eidors_obj('get-cache', inv_model, 'RtR_prior');
if ~isempty(RtR_prior)
   eidors_msg('calc_RtR_prior: using cached value', 3);
   return
end

if isfield(inv_model,'RtR_prior')
   RtR_prior= feval( inv_model.RtR_prior, inv_model );
elseif isfield(inv_model,'R_prior')
   % The user has provided an R prior. We can use this to
   % calculate RtR= R'*R;
   R= feval( inv_model.R_prior, inv_model );

   RtR_prior = R'*R;
else
   error('calc_RtR_prior: neither R_prior or RtR_prior func provided');
end

eidors_obj('set-cache', inv_model, 'RtR_prior', RtR_prior);
eidors_msg('calc_RtR_prior: setting cached value', 3);
