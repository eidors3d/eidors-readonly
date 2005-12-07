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

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_hyperparameter.m,v 1.1 2005-12-07 23:01:30 aadler Exp $

if isfield( inv_model.hyperparameter, 'func')

   hyperparameter = eidors_obj('get-cache', inv_model, 'hyperparameter');
   if ~isempty(hyperparameter)
      eidors_msg('calc_hyperparameter: using cached value', 2);
      return
   end

   hyperparameter= feval( inv_model.hyperparameter.func, inv_model);

   eidors_obj('set-cache', inv_model, 'hyperparameter', hyperparameter);
   eidors_msg('calc_hyperparameter: setting cached value', 2);

else
    hyperparameter= inv_model.hyperparameter.value;
end
