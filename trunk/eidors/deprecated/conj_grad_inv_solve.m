function img= conj_grad_inv_solve( varargin )
% CONJ_GRAD_INV_SOLVE inverse solver based on the CG
% inverse [Ref Shewchuck, 1994]
% img= conj_grad_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% parameters:
%   tol =     inv_model.parameters.term_tolerance;
%   maxiter = inv_model.parameters.max_iterations;

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','CONJ_GRAD_INV_SOLVE is deprecated as of 08-Jun-2012. Use INV_SOLVE_CONJ_GRAD instead.');
img = inv_solve_conj_grad(varargin{:});
