function Reg= noser_image_prior( inv_model );
% NOSER_IMAGE_PRIOR calculate image prior
% Reg= noser_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%
% Prior is diag( diag(J'*J)^exponent )
% param is normally .5, this value can be changed by
% setting inv_model.noser_image_prior.exponent= new_value

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','NOSER_IMAGE_PRIOR is deprecated as of 08-Jun-2012. Use PRIOR_NOSER instead.');
Reg = prior_noser(inv_model);
