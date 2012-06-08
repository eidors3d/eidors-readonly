function img= GN_abs_solve( varargin)
% GN_ABS_SOLVER absolute solver using Gauss Newton approximation
% img= gn_abs_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

% Step 1: fit to background
warning('EIDORS:deprecated','GN_ABS_SOLVE is deprecated as of 08-Jun-2012. Use INV_SOLVE_ABS_GN instead.');
img = inv_solve_abs_GN(varargin{:});
