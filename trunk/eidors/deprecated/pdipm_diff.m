function img=pdipm_diff( varargin )
% PDIPM_DIFF inverse solver for difference data using Primal/Dual interior point method
% img= ab_pdipm( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
%  inv_model.pdipm_diff.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.pdipm_diff.norm_image 1 or 2 (DEFAULT 2)
%  inv_model.pdipm_diff.beta     (default 1e-6)
%
% Parameters:
%  max_iters =  inv_model.parameters.max_iterations (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
% beta is the parameter that smooths the TV functional

% (C) 2008 Andrea Borsic. License: GPL version 2 or version 3
% $Id$


warning('EIDORS:deprecated','PDIPM_DIFF is deprecated as of 08-Jun-2012. Use INV_SOLVE_DIFF_PDIPM instead.');
img = inv_solve_diff_pdipm(varargin{:});
