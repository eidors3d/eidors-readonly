function Reg= np_calc_image_prior( inv_model );
% NP_CALC_IMAGE_PRIOR calculate image prior
% Ref= np_calc_image_prior( inv_model )
% Ref        => output regularization term
% inv_model  => inverse model struct

% $Id: np_calc_image_prior.m,v 1.1 2004-07-18 04:16:50 aadler Exp $

fwd_model= inv_model.fwd_model;

smooth_deg= inv_model.image_prior.parameters(1);
smooth_w  = inv_model.image_prior.parameters(2);

Reg = iso_f_smooth(fwd_model.elems, ...
                   fwd_model.nodes, ...
                   smooth_deg, smooth_w);

