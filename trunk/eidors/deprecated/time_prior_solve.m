function img= time_prior_solve( varargin )
% TIME_PRIOR_SOLVE inverse solver to account for time differences
% img= time_prior_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% TODO: This function really should be calling the proper
%   prior calculator functions, and not reimplementing
%   them internally

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','TIME_PRIOR_SOLVE is deprecated as of 08-Jun-2012. Use INV_SOLVE_TIME_PRIOR instead.');
img = inv_solve_time_prior(varargin{:});
