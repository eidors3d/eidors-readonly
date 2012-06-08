function img= backproj_solve( varargin)
% BACKPROJ_SOLVE inverse solver using backprojection
% NOTE: This is the beginnings of an attempt to reproduce
%  the backprojection algorithm implemented in the
%  Sheffield MKI EIT system. It is far from complete.
%
% If you wish to use the actual algorithm, use the
%  function "mk_common_gridmdl('backproj')"
%
% img= backproj_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% inv_model.backproj_solve.type = 'naive' (DEFAULT)
%    use naive (unfiltered algorithm)
% inv_model.backproj_solve.type = 'filtered' (NOT IMPLEMENTED YET)
%    ref: Barber DC Brown BH, "fast reconstruction of resistance
%         images", clin Phys Physiol Mes, pp 47-54, vol 8,sup A,1987
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','BACKPROJ_SOLVE is deprecated as of 08-Jun-2012. Use INV_SOLVE_BACKPROJ instead.');

img = inv_solve_backproj(varargin{:});
