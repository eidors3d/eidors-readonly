function RtR = calc_RtR_prior( inv_model, varargin )
% CALC_RtR_PRIOR: calculate image covariance R'*R
%   The image prior is matrix n_elem x n_elem of the
%   normalized a priori crosscorrelation FEM element values
% 
% calc_RtR_prior can be called as
%    RtR_prior= calc_RtR_prior( inv_model, ... )
%
% in each case it will call the inv_model.RtR_prior.func
%
% RtR_prior    the calculated RtR regularization prior
% inv_model    is an inv_model structure
%
% If a function to calculate RtR_prior is not provided,
% RtR = R_prior' * R_prior;

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_RtR_prior.m,v 1.2 2005-12-05 11:33:58 aadler Exp $

RtR_prior= eidors_obj('calc-or-cache', inv_model, ...
                 inv_model.RtR_prior.func, varargin{:} ); 

