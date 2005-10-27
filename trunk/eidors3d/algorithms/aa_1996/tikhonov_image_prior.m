function Reg= tikhonov_image_prior( inv_model );
% TIKHONOV_IMAGE_PRIOR calculate image prior
% Reg= tikhonov_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: tikhonov_image_prior.m,v 1.2 2005-10-27 13:28:08 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );

Reg = speye( pp.n_elem );

