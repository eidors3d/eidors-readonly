function R_prior = calc_R_prior( inv_model, varargin )
% R = calc_R_prior( inv_model, varargin )
% CALC_R_PRIOR: calculate regularization matrix R
%   The image prior is matrix n_elem x ??? 
% 
% calc_R_prior can be called as
%    R_prior= calc_R_prior( inv_model, ... )
%
% and will call the function inv_model.R_prior
% parameters to R_prior should be passed in the field
% inv_model.R_prior_function_name.parameters
% 
% If inv_model.R_prior is a matrix, calc_R_prior will return that matrix,
% possibly correcting for coarse2fine
%
% R_prior      calculated regularization prior R
% inv_model    is an inv_model structure

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id$


inv_model = rec_or_fwd_model( inv_model);


if isfield(inv_model,'R_prior')
   if isnumeric(inv_model.R_prior)
      R_prior = inv_model.R_prior;
   else
      R_prior= eidors_cache( inv_model.R_prior, inv_model );
   end
elseif isfield(inv_model,'RtR_prior')
   R_prior = eidors_cache(@calc_from_RtR_prior, inv_model,'calc_R_prior');
else
   error('calc_R_prior: neither R_prior or RtR_prior func provided');
end

if isfield(inv_model.fwd_model,'coarse2fine')
   c2f= inv_model.fwd_model.coarse2fine;
   if size(R_prior,1)==size(c2f,1)
%     we need to take into account coarse2fine - using a reasonable tol
      R_prior=R_prior*c2f;
   end
end

function R_prior = calc_from_RtR_prior(inv_model)
   % The user has provided an RtR prior. We can use this to
   % get R =RtR^(1/2). Not that this is non unique
   if isnumeric(inv_model.RtR_prior)
      RtR_prior = inv_model.RtR_prior;
   else
      RtR_prior= eidors_cache( inv_model.RtR_prior, inv_model );
   end
   
   % chol generates an error for rank deficient RtR_prior
   %     R_prior = chol (RtR_prior);
   % Instead we calculate cholinc with a droptol of 1e-5.
   %  For priors, this should be fine, since exact values
   %  especially far away, are not necessary
   ver = eidors_obj('interpreter_version');
   opts.droptol = 1e-5; 
   
   if ver.isoctave  || ver.ver < 7.012
      R_prior = cholinc(RtR_prior,opts.droptol);
   else
      R_prior = ichol(RtR_prior);
   end

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