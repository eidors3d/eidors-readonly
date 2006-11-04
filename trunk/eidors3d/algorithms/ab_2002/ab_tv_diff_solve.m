function img= ab_tv_diff_solve( inv_model, data1, data2)
% AB_TV_DIFF_SOLVE inverse solver for Andrea Borsic's
%   Total Variation solver for use with difference EIT
% img= ab_tv_diff_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: ab_tv_diff_solve.m,v 1.6 2006-11-04 15:07:39 aadler Exp $


[alpha1,alpha2,beta,maxiter,tol,keepiters]= get_params(inv_model);

try
   normalize = inv_model.fwd_model.normalize_measurements;
catch
   normalize = 0;
end

if normalize
   dva= 1 - data2 ./ data1;
else   
   dva= data1 - data2;
end


sol=primaldual_tvrecon_lsearch(inv_model, dva ,maxiter,alpha1,alpha2);
if ~keepiters
   sol=sol(:,end);
end

img.name= 'solved by ab_tv_diff_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;

function [alpha1,alpha2,beta,maxiter,tol,keepiters]= ...
          get_params(inv_model);
   alpha1= calc_hyperparameter( inv_model);
   alpha2= 1e-4;
   beta= 1;
   if isfield(inv_model,'ab_tv_diff_solve')
      if isfield(inv_model.ab_tv_diff_solve,'alpha2')
         alpha2= inv_model.ab_tv_diff_solve.alpha2;
      end
      if isfield(inv_model.ab_tv_diff_solve,'beta')
         alpha2= inv_model.ab_tv_diff_solve.beta;
      end
   end

   maxiter= 3;
   tol= 1e-4;
   keepiters= 0;
   if isfield(inv_model,'parameters')
      tol =     inv_model.parameters.term_tolerance;
      maxiter = inv_model.parameters.max_iterations;

      if isfield(inv_model.parameters,'keep_iterations')
         keepiters = inv_model.parameters.keep_iterations;
      end
   end
