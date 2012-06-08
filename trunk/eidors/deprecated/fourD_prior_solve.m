function img= fourD_prior_solve( varargin )
% fourD_prior_solve-- inverse solver to account for temporal
% and 3D spatial correlation
% img= fourD_prior_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2007, Tao Dai and Andy Adler. Licenced under the GPL Version 2
% $Id$

warning('EIDORS:deprecated','FOURD_PRIOR_SOLVE is deprecated as of 08-Jun-2012. Use INV_SOLVE_4D_PRIOR instead.');
img = inv_solve_4d_prior(varargin{:});
