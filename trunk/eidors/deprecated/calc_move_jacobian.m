function J = calc_move_jacobian(varargin)
% CALC_MOVE_JACOBIAN   Computes the Jacobian matrix for conductivity and
% electrode movement variables in 3D EIT.
% Args:     fwd_model - the EIDORS object forward model
%            img_bkgd - the image background conductivity
% Returns:          J - the Jacobian matrix [Jc, Jm]
%
% WARNING: THIS CODE IS EXPERIMENTAL AND GIVES PROBLEMS
% SEE: Camille Gomez-Laberge, Andy Adler
% Direct EIT Jacobian calculations for conductivity change
%  and electrode movement,  Physiol. Meas., 29:S89-S99, 2008

% (C) 2007, Camille Gomez-Laberge and Andy Adler.
%  License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','CALC_MOVE_JACOBIAN is deprecated as of 08-Jun-2012. Use JACOBIAN_MOVEMENT instead.');

J = jacobian_movement(varargin{:});
