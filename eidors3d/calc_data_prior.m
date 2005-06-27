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
% Many common EIT (and other regularized) algorithms do not
%  contain a data prior term. For these algorithms, this function
%  generates a reasonable approximation based on uniform noise.
%
% $Id: calc_data_prior.m,v 1.3 2005-06-27 15:55:51 aadler Exp $

data_prior = eidors_obj('get-cache', inv_model, 'data_prior');

if ~isempty(data_prior)
   eidors_msg('calc_data_prior: using cached value', 2);
   return
end

% work around the change in meaning of && in matlab 6 - 6.5
data_prior= [];
if isfield(inv_model,'data_prior')
if isfield(inv_model.data_prior,'func')
if ~isempty('inv_model.data_prior.func')
      data_prior= feval( inv_model.data_prior.func, inv_model);
end
end
end

if isempty(data_prior)
   data_prior= default_data_prior( inv_model );
end

eidors_obj('set-cache', inv_model, 'data_prior', data_prior);
eidors_msg('calc_data_prior: setting cached value', 2);
 
% Calculate a data prior for an assumption of uniform noise
% on each channel
% 
function data_prior = default_data_prior( inv_model )

   fwd_model= inv_model.fwd_model;
   if isfield(fwd_model,'normalize_measurements')
      normalize = fwd_model.normalize_measurements;
   else
      normalize = 0;
   end

   n =  calc_n_meas( fwd_model );

   if ~normalize
      data_prior= speye( n );
   else
      % How to do this??
      data_prior= speye( n );
   end

function n_meas = calc_n_meas( fwd_model )

   n_meas = 0;
   for i= 1:length(fwd_model.stimulation );
       n_meas = size(fwd_model.stimulation(i).meas_pattern,??);
   end
