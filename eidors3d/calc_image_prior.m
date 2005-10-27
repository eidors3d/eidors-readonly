function image_prior = calc_image_prior( inv_model, varargin )
% CALC_DATA_PRIOR: calculate prior probabilities for image term
%   The image prior is matrix n_elem x n_elem of the a priori
%     crosscorrelation FEM element values
% 
% calc_image_prior can be called as
%    image_prior= calc_image_prior( inv_model, ... )
%
% in each case it will call the inv_model.image_prior.func
%
% image_prior   is the calculated image prior
% inv_model    is an inv_model structure

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_image_prior.m,v 1.8 2005-10-27 13:28:08 aadler Exp $

image_prior= eidors_obj('calc-or-cache', inv_model, ...
                 inv_model.image_prior.func, varargin{:} ); 

