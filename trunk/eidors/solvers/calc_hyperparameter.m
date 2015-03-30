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

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isfield( inv_model.hyperparameter, 'func')
   try inv_model.hyperparameter = str2func(inv_model.hyperparameter); end
   hyperparameter = eidors_cache(...
                    inv_model.hyperparameter.func,{inv_model},'hyperparameter');
else
    hyperparameter= inv_model.hyperparameter.value;
end
