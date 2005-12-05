function RtR_prior = calc_RtR_prior( inv_model, varargin )
% RtR = calc_RtR_prior( inv_model, varargin )
% CALC_RtR_PRIOR: calculate image covariance R'*R
%   The image prior is matrix n_elem x n_elem of the
%   normalized a priori crosscorrelation FEM element values
% 
% calc_RtR_prior can be called as
%    RtR_prior= calc_RtR_prior( inv_model, ... )
%
% in each case it will call the inv_model.RtR_prior.func
%
% RtR_prior    the calculated RtR regularization prior
% inv_model    is an inv_model structure
%
% If a function to calculate RtR_prior is not provided,
% RtR = R_prior' * R_prior;

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_RtR_prior.m,v 1.4 2005-12-05 22:12:11 aadler Exp $

RtR_prior = eidors_obj('get-cache', inv_model, 'RtR_prior');
if ~isempty(RtR_prior)
   eidors_msg('calc_RtR_prior: using cached value', 2);
   return
end

if isfield(inv_model,'RtR_prior')
   RtR_prior= feval( inv_model.RtR_prior.func, inv_model );
elseif isfield(inv_model,'R_prior')
   % The user has provided an R prior. We can use this to
   % calculate RtR= R'*R;
   R= feval( inv_model.R_prior.func, inv_model );

   RtR_prior = R'*R;
else
   error('calc_RtR_prior: neither R_prior or RtR_prior func provided');
end

eidors_obj('set-cache', inv_model, 'RtR_prior', RtR_prior);
eidors_msg('calc_RtR_prior: setting cached value', 2);



