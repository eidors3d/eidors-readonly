function Reg= tikhonov_image_prior( inv_model );
% TIKHONOV_IMAGE_PRIOR calculate image prior
% Reg= tikhonov_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

pp= aa_fwd_parameters( inv_model.fwd_model );

Reg = speye( pp.n_elem );

