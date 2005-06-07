function Reg= tikhonov_image_prior( inv_model );
% TIKHONOV_IMAGE_PRIOR calculate image prior
% Reg= tikhonov_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% $Id: tikhonov_image_prior.m,v 1.1 2005-06-07 02:47:31 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );

Reg = speye( pp.n_elem );

