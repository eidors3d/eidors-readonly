function ok= calc_model_prior_test;
% Verify model prior calcs
% $Id$

% TODO: also test the various inverse prior calls

imdl= mk_common_model('c2c2',16);
try; imdl= rmfield(imdl,'RtR_prior'); end
try; imdl= rmfield(imdl,'R_prior');   end

any_priors= {@prior_tikhonov, ...
             @prior_noser, ...
             @prior_gaussian_HPF, ...
             @prior_laplace};

R_priors=   {any_priors{:}, ...
             @prior_TV};

% Call R_priors as R_priors
eidors_cache clear
for i = 1:length(R_priors); p = R_priors{i};
   inv_mdl= imdl;
   inv_mdl.R_prior= p;
   R= calc_R_prior(inv_mdl);
   fprintf('R_prior: %20s  R_condest= %5.4g\n', func2str(p), condest(R'*R));
end
   
% Call R_priors as RtR_priors
eidors_cache clear
for i = 1:length(R_priors); p = R_priors{i};
   inv_mdl= imdl;
   inv_mdl.R_prior= p;
   RtR= calc_RtR_prior(inv_mdl);
   fprintf('R_prior: %20s  RtR_condest= %5.4g\n', func2str(p), condest(RtR));
end
   
% Call RtR_priors as RtR_priors
eidors_cache clear
for i = 1:length(R_priors); p = R_priors{i};
   inv_mdl= imdl;
   inv_mdl.RtR_prior= p;
   RtR= calc_RtR_prior(inv_mdl);
   if diff(size(RtR))~=0  % non-square
      fprintf('RtR_prior: %20s  RtR_condest= NON-SQUARE\n', func2str(p) );
   else
      fprintf('RtR_prior: %20s  RtR_condest= %5.4g\n', func2str(p), condest(RtR));
   end
end

% Call RtR_priors as R_priors
eidors_cache clear
for i = 1:length(R_priors); p = R_priors{i};
   inv_mdl= imdl;
   inv_mdl.RtR_prior= p;
   if strcmp(func2str(p), 'prior_TV')
      continue; % not fair to ask it to ichol a non-square matrix
   end
   R= calc_R_prior(inv_mdl);
   fprintf('RtR_prior: %20s  R_condest= %5.4g\n', func2str(p), condest(R'*R));
end
