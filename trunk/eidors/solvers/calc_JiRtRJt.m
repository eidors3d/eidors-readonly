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
% if JiRtRJt is a matrix, than calc_JiRtRJt will return that matrix
%
% JiRtRJt         calculated regularization prior R
% inv_model       is an inv_model structure
% inv_model.JiRtRJt_func function to make calculation
%
% TODO: think about how to implement this!!

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isfield(inv_model,'JiRtRJt') 
   if isnumeric(inv_model.JiRtRJt)
      JiRtRJt = inv_model.JiRtRJt;
   else
      try inv_model.JiRtRJt = str2func(inv_model.JiRtRJt); end
      JiRtRJt = eidors_cache(inv_model.JiRtRJt,{inv_model});
   end
   return
end

% try to calculate from scratch
JiRtRJt = eidors_cache(@calc_JiRtRJt_from_scratch,{inv_model},'calc_JiRtRJt');


function JiRtRJt = calc_JiRtRJt_from_scratch(inv_model)
   eidors_msg( ...
      'calc_JiRtRJt: trying to calculate JiRtRJt from RtR_prior',2);
   RtR_prior = calc_RtR_prior(inv_model);
   % regularize slightly so inverse exists
   % This is very slow, but I can't think of a better idea
   RtR_p_reg = spdiags( spdiags(RtR_prior,0)*1.00001, 0, RtR_prior);

   img_bkgnd= calc_jacobian_bkgnd( inv_model );
   J = calc_jacobian( img_bkgnd);

   JiRtRJt= J*(RtR_p_reg\J');
