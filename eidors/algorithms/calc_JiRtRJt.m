function [JiRtRJt,iRtRJt] = calc_JiRtRJt( inv_model, varargin )
% [JiRtRJt,iRtRJt] = calc_JiRtRJt( inv_model, varargin )
% CALC_iRtR_PRIOR: calculate regularization matrix J*inv(R'*R)*J'
%   This is a model of the covariance of image elements
%   The image prior is matrix n_elem x n_elem 
% 
% calc_JiRtRJt can be called as
%    JiRtRJt= calc_JiRtRJt( inv_model, ... )
%
% and will call the function inv_model.JiRtRJt
% parameters to JiRtRJt should be passed in the field
% inv_model.JiRtRJt_function_name.parameters
%
% JiRtRJt         calculated regularization prior R
% inv_model       is an inv_model structure
% inv_model.JiRtRJt_func function to make calculation
%
% TODO: think about how to implement this!!

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

JiRtRJt = eidors_obj('get-cache', inv_model, 'JiRtRJt');
if ~isempty(JiRtRJt)
   eidors_msg('calc_JiRtRJt: using cached value', 3);
   return
end

if isfield(inv_model,'JiRtRJt')
   JiRtRJt= feval( inv_model.JiRtRJt, inv_model );
else
   eidors_msg( ...
      'calc_JiRtRJt: trying to calculate JiRtRJt from RtR_prior',2);
   RtR_prior = feval( inv_model.RtR_prior, inv_model );
   % regularize slightly so inverse exists
   % This is very slow, but I can't think of a better idea
   RtR_p_reg = spdiags( spdiags(RtR_prior,0)*1.00001, 0, RtR_prior);

   img_bkgnd= calc_jacobian_bkgnd( inv_model );
   J = calc_jacobian( fwd_model, img_bkgnd);

   JiRtRJt= J*(RtR_p_reg\J');
end

eidors_obj('set-cache', inv_model, 'JiRtRJt', JiRtRJt);
eidors_msg('calc_JiRtRJt: setting cached value', 3);
