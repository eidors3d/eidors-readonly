function Reg= np_calc_image_prior( inv_model );
% NP_CALC_IMAGE_PRIOR calculate image prior
% Ref= np_calc_image_prior( inv_model )
% Ref        => output regularization term
% inv_model  => inverse model struct

% $Id: np_calc_image_prior.m,v 1.2 2005-02-23 14:43:46 aadler Exp $

Reg = eidors_obj('cache', inv_model, 'image_prior');
if ~isempty(Reg)
   eidors_msg('np_calc_image_prior: using cached value',1);
   return
end

smooth_deg= inv_model.image_prior.parameters(1);
smooth_w  = inv_model.image_prior.parameters(2);

Reg = iso_f_smooth(inv_model.fwd_model.elems, ...
                   inv_model.fwd_model.nodes, ...
                   smooth_deg, smooth_w);

eidors_obj('cache', inv_model, 'image_prior', Reg);
eidors_msg('np_calc_image_prior: setting cached value',1);
