function J= aa_calc_jacobian( varargin );
% AA_CALC_JACOBIAN: J= aa_calc_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for current stimulation EIT
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','AA_CALC_JACOBIAN is deprecated as of 06-Jun-2012. Use CALC_JACOBIAN_ADJOINT instead.');

J = calc_jacobian_adjoint( varargin{:} );

