function img= aa_inv_conj_grad( varargin )
% AA_INV_CONJ_GRAD inverse solver based on the CG
% inverse [Ref Shewchuck, 1994]
% img= aa_inv_conj_grad( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','AA_INV_CONJ_GRAD is deprecated as of 07-Jun-2012. Use CONJ_GRAD_INV_SOLVE instead.');

img = conj_grad_inv_solve( varargin{:} );

