function img= inv_kalman_diff( inv_model, data1, data2)
% INV_KALMAN_DIFF inverse solver for difference EIT
% img= inv_kalman_diff( inv_model, data1, data2)
%
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% if inv_model.fwd_model.stimulation(:).delta_time
%   exists and is non_zero, then the kalman filter will
%   be applied to each data measurement separately
%
% Note that the classic Kalman filter assumes that the
%   time step between each measurement is constant
%   (ie as part of the state update eqn). inv_kalman_diff
%   cannot work with non-constant time steps
%
% if inv_model.inv_kalman_diff.keep_K_k1 = 1
%  then img outputs img.inv_kalman_diff.K_k1 = K_k1
%  this can be used to estimate noise properties

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','INV_KALMAN_DIFF is deprecated as of 08-Jun-2012. Use INV_SOLVE_DIFF_KALMAN instead.');
img = inv_solve_diff_kalman(varargin{:});
