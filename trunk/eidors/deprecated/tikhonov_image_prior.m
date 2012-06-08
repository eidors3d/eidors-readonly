function Reg= tikhonov_image_prior( inv_model );
% TIKHONOV_IMAGE_PRIOR calculate image prior
% Reg= tikhonov_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

%pp= fwd_model_parameters( inv_model.fwd_model );

warning('EIDORS:deprecated','TIKHONOV_IMAGE_PRIOR is deprecated as of 08-Jun-2012. Use PRIOR_TIKHONOV instead.');
Reg = prior_tikhonov(inv_model);
