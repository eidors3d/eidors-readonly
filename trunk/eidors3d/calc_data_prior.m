function data_prior = calc_data_prior( inv_model )
% CALC_DATA_PRIOR: calculate prior probabilities for data term
%   The data prior is matrix n_meas x n_meas of the a priori
%     crosscorrelation of measurements
% 
% calc_data_prior can be called as
%    data_prior= calc_data_prior( inv_model )
%
% in each case it will call the inv_model.data_prior.func
%
% data_prior   is the calculated data prior
% inv_model    is an inv_model structure
%
% $Id: calc_data_prior.m,v 1.2 2005-01-20 21:26:29 billlion Exp $

data_prior= feval( inv_model.data_prior.func, inv_model);
 