function Reg= aa_e_move_image_prior( inv_model );
% AA_E_MOVE_IMAGE_PRIOR calculate image prior
% Reg= aa_e_move_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   inv_model.image_prior.parameters(1) -> relative weighting
%     of movement vs image fraction of hyperparameter
%     => Default = 100
%   inv_model.aa_e_move_image_prior.RegC.func = Cond Reg fcn
%   inv_model.aa_e_move_image_prior.RegM.func = Move Reg fcn
%   either @laplace_movement_image_prior OR @tikhonov_movement_image_prior
%
% For image portion, we use a Laplace prior, as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself
%
% For the movmenent portion, we define a smoothness
% constraint, such that Rij = -1 for adjacent electrodes
%
% If used with a dual model (ie coarse2fine mapping), ensure
%    imdl.prior_use_fwd_not_rec = 1;

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% relative strengths of conductivity and movement priors
warning('EIDORS:deprecated','AA_E_MOVE_IMAGE_PRIOR is deprecated as of 07-Jun-2012. Use PRIOR_MOVEMENT instead.');

Reg = prior_movement( inv_model);
