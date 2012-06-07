function J= aa_e_move_jacobian( varargin );
% AA_E_MOVE_JACOBIAN: J= aa_e_move_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for EIT, based on conductivity
%   change and movement of electrodes
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.conductivity_jacobian = fcn
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','AA_E_MOVE_JACOBIAN is deprecated as of 07-Jun-2012. Use CALC_MOVE_JACOBIAN_PERTURB instead.');

J = calc_move_jacobian_perturb( varargin{:} );
