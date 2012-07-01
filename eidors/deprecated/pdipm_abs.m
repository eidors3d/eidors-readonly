function img=pdipm_abs( varargin )
% PDIPM_ABS  inverse solver for absolute data using Primal/Dual interior point method
% img= pdipm_abs( inv_model, data);
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data       => vector of eit data
%
%  inv_model.pdipm_abs.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.pdipm_abs.norm_prior 1 or 2 (DEFAULT 2)
%  inv_model.pdipm_abs.beta     (default 1e-6)
%
% Parameters:
%  max_iter =  inv_model.parameters.max_iteration (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
% beta is the parameter that smooths the TV functional

% (C) 2010 Andrea Borsic + Andy Adler. License: GPL v2 or v3
% $Id$


warning('EIDORS:deprecated','PDIPM_ABS is deprecated as of 08-Jun-2012. Use INV_SOLVE_ABS_PDIPM instead.');

if isfield(inv_model,'pdipm_abs');
  inv_model.inv_solve_abs_pdipm = inv_model.pdipm_abs;
end

img = inv_solve_abs_pdipm(varargin{:});
