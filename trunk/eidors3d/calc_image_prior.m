function image_prior = calc_image_prior( inv_model )
% CALC_DATA_PRIOR: calculate prior probabilities for image term
%   The image prior is matrix n_elem x n_elem of the a priori
%     crosscorrelation FEM element values
% 
% calc_image_prior can be called as
%    image_prior= calc_image_prior( inv_model )
%
% in each case it will call the inv_model.image_prior.func
%
% image_prior   is the calculated image prior
% inv_model    is an inv_model structure
%
% $Id: calc_image_prior.m,v 1.1 2004-07-18 04:16:51 aadler Exp $

image_prior= feval( inv_model.image_prior.func, inv_model);
