function Reg= gaussian_HPF_prior( inv_model );
% GAUSSIAN_HPF_PRIOR calculate image prior
% Reg= gaussian_HPF_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   diam_frac= inv_model.fwd_model.gaussian_HPF_prior.diam_frac DEFAULT 0.1

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','GAUSSIAN_HPF_PRIOR is deprecated as of 08-Jun-2012. Use PRIOR_GAUSSIAN_HPF instead.');
prior_gaussian_hpf(inv_model);
