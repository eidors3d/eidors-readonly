function R_prior = calc_R_prior( inv_model, varargin )
% R = calc_R_prior( inv_model, varargin )
% CALC_R_PRIOR: calculate regularization matrix R
%   The image prior is matrix n_elem x ??? 
% 
% calc_R_prior can be called as
%    R_prior= calc_R_prior( inv_model, ... )
%
% in each case it will call the inv_model.image_prior.func
%
% R_prior      calculated regularization prior R
% inv_model    is an inv_model structure

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_R_prior.m,v 1.3 2005-12-05 13:13:31 aadler Exp $

R_prior= eidors_obj('calc-or-cache', inv_model, ...
                 inv_model.R_prior.func, varargin{:} ); 

