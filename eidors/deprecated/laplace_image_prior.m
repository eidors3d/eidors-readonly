function Reg= laplace_image_prior( inv_model );
% LAPLACE_IMAGE_PRIOR calculate image prior
% Reg= laplace_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%
% This image prior is intended to be used as
%  R'*R, but may be used as R for as well.
%
% The Laplace prior is a 2nd order high pass filter.
% On a rectangular mesh, it is a convolution with
%   [-1,-1,-1;      [ 0;-1; 0
%    -1, 8,-1    or  -1; 4;-1
%    -1,-1,-1]        0;-1; 0]
%
% On a finite element mesh, we define the it as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','LAPLACE_IMAGE_PRIOR is deprecated as of 08-Jun-2012. Use PRIOR_LAPLACE instead.');
Reg = prior_laplace(inv_model);
