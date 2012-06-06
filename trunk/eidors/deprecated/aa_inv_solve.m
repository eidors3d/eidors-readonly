function img= GN_one_step_diff_solve( varargin )
% AA_INV_SOLVE inverse solver using approach of Adler&Guardo 1996
% img= GN_one_step_diff_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','AA_INV_SOLVE is deprecated as of 06-Jun-2012. Use GN_ONE_STEP_DIFF_SOLVE instead.');

img = GN_one_step_diff_solve( varargin{:} );
