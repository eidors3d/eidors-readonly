function ok= calc_model_prior_test;
% Verify model prior calcs
% $Id: calc_model_prior_test.m,v 1.1 2008-03-18 20:44:16 aadler Exp $

% TODO: also test the various inverse prior calls

imdl= mk_common_model('c2c2',16);
try; imdl= rmfield(imdl,'RtR_prior'); end
try; imdl= rmfield(imdl,'R_prior');   end

any_priors= {@tikhonov_image_prior, ...
             @noser_image_prior, ...
             @Gaussian_HPF_prior, ...
             @laplace_image_prior};

R_priors=   {any_priors{:}, ...
             @ab_calc_tv_prior};

% Call R_priors as R_priors
eidors_cache clear
for p = [R_priors{:}]
   inv_mdl= imdl;
   inv_mdl.R_prior= p;
   R= calc_R_prior(inv_mdl);
   fprintf('R_prior: %20s  R_condest= %5.4g\n', func2str(p), condest(R'*R));
end
   
% Call R_priors as RtR_priors
eidors_cache clear
for p = [R_priors{:}]
   inv_mdl= imdl;
   inv_mdl.R_prior= p;
   RtR= calc_RtR_prior(inv_mdl);
   fprintf('R_prior: %20s  RtR_condest= %5.4g\n', func2str(p), condest(RtR));
end
   
% Call RtR_priors as RtR_priors
eidors_cache clear
for p = [any_priors{:}]
   inv_mdl= imdl;
   inv_mdl.RtR_prior= p;
   RtR= calc_RtR_prior(inv_mdl);
   fprintf('RtR_prior: %20s  RtR_condest= %5.4g\n', func2str(p), condest(RtR));
end

% Call RtR_priors as R_priors
eidors_cache clear
for p = [any_priors{:}]
   inv_mdl= imdl;
   inv_mdl.RtR_prior= p;
   R= calc_R_prior(inv_mdl);
   fprintf('RtR_prior: %20s  R_condest= %5.4g\n', func2str(p), condest(R'*R));
end
