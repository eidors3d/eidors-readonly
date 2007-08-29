function Reg= np_calc_image_prior( inv_model );
% NP_CALC_IMAGE_PRIOR calculate image prior
%
% Intended to be used as calc_R_prior
%
% Ref= np_calc_image_prior( inv_model )
% Ref        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: np_calc_image_prior.m,v 1.5 2007-08-29 09:00:56 aadler Exp $

parameters=  inv_model.np_calc_image_prior.parameters;
smooth_deg= parameters(1);
smooth_w  = parameters(2);

Reg = iso_f_smooth(inv_model.fwd_model.elems, ...
                   inv_model.fwd_model.nodes, ...
                   smooth_deg, smooth_w);
