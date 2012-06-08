function Reg= calc_covar_prior( inv_model )
% CALC_COVAR_PRIOR image prior with distance-based interelement covar
% This is a simplification of exponential_covar_prior.m
% Reg= calc_covar_prior( inv_model )
% Reg        => output regularization 
% inv_model  => inverse model struct
% P_type--prior type
% 1: elements are globally correlated
% 2: elements within/without electrode rings are correlated to elements in same region.
% 3: only elements within electrode rings are correlated.

% (C) 2007, Tao Dai and Andy Adler. Licenced under the GPL Version 2
% $Id$

% get average x,y,z of each element
warning('EIDORS:deprecated','CALC_COVAR_PRIOR is deprecated as of 08-Jun-2012. Use PRIOR_COVAR instead.');

Reg = prior_covar(inv_model);
