function likelihood= gaussian_prior_likelyhood( varargin )
% Parameters for image
%   inv_model.gaussian_prior_likelihood.img_mean -> image mean
%   inv_model.gaussian_prior_likelihood.R_prior -> L*L' = inv(image covariance)
%   inv_model.gaussian_prior_likelihood.img_exp -> ( default = 2)
% Parameters for data
%   inv_model.gaussian_prior_likelihood.Noise -> L*L' = inv(Noise covariance)
%   inv_model.gaussian_prior_likelihood.data_exp -> ( default = 2)
%
% Function to be used with mcmc_solve

% (C) 2007 Nick Polydorides. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','GAUSSIAN_PRIOR_LIKELYHOOD is deprecated as of 08-Jun-2012. Use PRIOR_GAUSSIAN_LIKELIHOOD instead.');

if isfield(inv_model,'gaussian_prior_likelyhood');
  inv_model.prior_gaussian_likelihood = inv_model.gaussian_prior_likelyhood;
end

likelihood = prior_gaussian_likelihood(varargin{:});
