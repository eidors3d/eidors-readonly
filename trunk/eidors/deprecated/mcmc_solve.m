function img= mcmc_solve( varargin )
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
% inv_model.parameters.max_iterations (default 1);
% inv_model.parameters.term_tolerance (default 1e-3);

% (C) 2007 Nick Polydorides. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','MCMC_SOLVE is deprecated as of 08-Jun-2012. Use INV_SOLVE_MCMC instead.');
img = inv_solve_mcmc(varargin{:});
