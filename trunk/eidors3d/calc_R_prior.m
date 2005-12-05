function R_prior = calc_R_prior( inv_model, varargin )
% R = calc_R_prior( inv_model, varargin )
% CALC_R_PRIOR: calculate regularization matrix R
%   The image prior is matrix n_elem x ??? 
% 
% calc_R_prior can be called as
%    R_prior= calc_R_prior( inv_model, ... )
%
% in each case it will call the inv_model.R_prior.func
%
% R_prior      calculated regularization prior R
% inv_model    is an inv_model structure

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_R_prior.m,v 1.4 2005-12-05 22:12:11 aadler Exp $

R_prior = eidors_obj('get-cache', inv_model, 'R_prior');
if ~isempty(R_prior)
   eidors_msg('calc_R_prior: using cached value', 2);
   return
end

if isfield(inv_model,'R_prior')
   R_prior= feval( inv_model.R_prior.func, inv_model );
elseif isfield(inv_model,'RtR_prior')
   % The user has provided an RtR prior. We can use this to
   % get R =RtR^(1/2). Not that this is non unique
   RtR_prior= feval( inv_model.RtR_prior.func, inv_model );

   % generates an error for rank deficient RtR_prior
   R_prior = chol (RtR_prior);
else
   error('calc_R_prior: neither R_prior or RtR_prior func provided');
end

eidors_obj('set-cache', inv_model, 'R_prior', R_prior);
eidors_msg('calc_R_prior: setting cached value', 2);



