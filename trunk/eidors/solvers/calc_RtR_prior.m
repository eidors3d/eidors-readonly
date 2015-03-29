function RtR_prior = calc_RtR_prior( inv_model )
% RtR = calc_RtR_prior( inv_model )
% CALC_RtR_PRIOR: calculate image regularization prior
%   R'*R (which is an estimate of the inverse of the covariance)
%
%   Typically, the image prior is matrix n_elem x n_elem of the
%   normalized a priori crosscorrelation FEM element values
% 
% calc_RtR_prior can be called as
%    RtR_prior= calc_RtR_prior( inv_model, ... )
%
% and will call the function inv_model.RtR_prior
% parameters to RtR_prior should be passed in the field
% inv_model.RtR_prior_function_name.parameters
%
% If inv_model.RtR_prior is a matrix, calc_RtR_prior will return that matrix,
% possibly correcting for coarse2fine
%
% if there exists a field inv_model.rec_model, then
%   the prior is calculated on the rec_model rather than
%   the fwd_model. This will not be done if 
% inv_model.prior_use_fwd_not_rec= 1;
%
% RtR_prior    the calculated RtR regularization prior
% inv_model    is an inv_model structure
%
% If a function to calculate RtR_prior is not provided,
% RtR = R_prior' * R_prior;

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

inv_model = rec_or_fwd_model( inv_model);

if isfield(inv_model,'RtR_prior')
   if isnumeric(inv_model.RtR_prior)
      RtR_prior = inv_model.RtR_prior;
   else
      RtR_prior= eidors_cache( inv_model.RtR_prior, inv_model );
   end
elseif isfield(inv_model,'R_prior')
   RtR_prior = eidors_cache(@calc_from_R_prior, inv_model, 'calc_RtR_prior');
else
   error('calc_RtR_prior: neither R_prior nor RtR_prior provided');
end

if isfield(inv_model.fwd_model,'coarse2fine')
   c2f= inv_model.fwd_model.coarse2fine;
   if size(RtR_prior,1)==size(c2f,1)
%     we need to take into account coarse2fine - using a reasonable tol
      eidors_msg('calc_RtR_prior: using coarse2fine to model RtR');
      f2c= c2f'; %pinv(c2f,1e-6);
      RtR_prior=f2c*RtR_prior*c2f;
   end
end

function RtR_prior = calc_from_R_prior(inv_model)

   % The user has provided an R prior. We can use this to
   % calculate RtR= R'*R;
   if isnumeric(inv_model.R_prior)
      R = inv_model.R_prior;
   else
      R= eidors_cache( inv_model.R_prior, inv_model );
   end

   RtR_prior = R'*R;
   
   

function inv_model = rec_or_fwd_model( inv_model);

   if isfield(inv_model,'rec_model');
      use_rec_model = 1;
      try if inv_model.prior_use_fwd_not_rec== 1;
         use_rec_model = 0;
      end; end

      if use_rec_model
          % copy the normalize flag from the fwd_model to prevent warnings
         inv_model.rec_model = mdl_normalize(inv_model.rec_model, ...
                                       mdl_normalize(inv_model.fwd_model));
         inv_model.fwd_model= inv_model.rec_model;
         inv_model= rmfield(inv_model,'rec_model');
      end
   end
