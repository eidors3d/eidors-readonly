function hyperparameter = calc_hyperparameter( inv_model )
% CALC_HYPERPARAMETER: calculate hyperparameter value
%   The hyperparameter is can be either provided directly,
%     or can be based an automatic selection approach
% 
% calc_hyperparameter can be called as
%    hyperparameter= calc_hyperparameter( inv_model )
% where inv_model    is an inv_model structure
%
% if inv_model.hyperparameter.func exists, it will be
%   called, otherwise inv_model.hyperparameter.value will
%   be returned
%
% TODO: does hyperparameter depend on inv_model, or does
%       it also depend on the data?
%
% $Id: calc_hyperparameter.m,v 1.1 2005-06-07 03:06:58 aadler Exp $

hyperparameter = eidors_obj('get-cache', inv_model, 'hyperparameter');

if ~isempty(hyperparameter)
   eidors_msg('calc_image_prior: using cached value', 2);
   return
end

if isfield( inv_model.hyperparameter, 'func')
    hyperparameter= feval( inv_model.hyperparameter.func, inv_model);
else
    hyperparameter= inv_model.hyperparameter.value;
end

eidors_obj('set-cache', inv_model, 'hyperparameter', hyperparameter);
eidors_msg('calc_hyperparameter: setting cached value', 2);
