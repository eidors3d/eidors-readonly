function J= perturb_jacobian( varargin )
% PERTURB_JACOBIAN: J= perturb_jacobian( fwd_model, img)
% Calculate Jacobian Matrix, based on small perturbations
%   in the forward model. This will tend to be slow, but
%   should be best used to 'sanity check' other code
%
% J         = Jacobian matrix
% fwd_model = forward model
% fwd_model.perturb_jacobian.delta   - delta perturbation to use
% img = image background for jacobian calc

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','PERTURB_JACOBIAN is deprecated as of 08-Jun-2012. Use JACOBIAN_PERTURB instead.');

if isfield(inv_model,'perturb_jacobian');
  inv_model.jacobian_perturb = inv_model.perturb_jacobian;
end

J = jacobian_perturb(varargin{:});
