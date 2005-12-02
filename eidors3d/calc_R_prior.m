function image_prior = calc_image_prior( inv_model, varargin )
% CALC_R_PRIOR: calculate regularization matrix R
%   The image prior is matrix n_elem x ??? 
% 
% calc_image_prior can be called as
%    image_prior= calc_image_prior( inv_model, ... )
%
% in each case it will call the inv_model.image_prior.func
%
% image_prior   is the calculated image prior
% inv_model    is an inv_model structure

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_R_prior.m,v 1.1 2005-12-02 11:59:17 aadler Exp $

image_prior= eidors_obj('calc-or-cache', inv_model, ...
                 inv_model.image_prior.func, varargin{:} ); 

