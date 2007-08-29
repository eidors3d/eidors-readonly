function Reg= tikhonov_image_prior( inv_model );
% TIKHONOV_IMAGE_PRIOR calculate image prior
% Reg= tikhonov_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: tikhonov_image_prior.m,v 1.8 2007-08-29 09:19:24 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );

Reg = speye( pp.n_elem );

