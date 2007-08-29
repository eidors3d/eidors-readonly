function img= ab_tv_diff_solve( inv_model, data1, data2)
% AB_TV_DIFF_SOLVE inverse solver for Andrea Borsic's
%   Total Variation solver for use with difference EIT
% img= ab_tv_diff_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: ab_tv_diff_solve.m,v 1.5 2007-08-29 09:15:30 aadler Exp $


[alpha1,alpha2,beta,maxiter,tol,keepiters]= get_params(inv_model);

if inv_model.fwd_model.normalize_measurements;
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end


for i=1:size(dva,2)
   soln=primaldual_tvrecon_lsearch(inv_model, dva(:,i), ...
        maxiter,alpha1,alpha2, tol, beta);
   if ~keepiters & size(dva,2)>1
      sol(:,i)=soln(:,end);
   else
      sol=soln;
   end
end

img.name= 'solved by ab_tv_diff_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;

function [alpha1,alpha2,beta,maxiter,tol,keepiters]= ...
          get_params(inv_model);
   alpha1= 1e-2;
   beta= .01;
   if isfield(inv_model,'ab_tv_diff_solve')
      if isfield(inv_model.ab_tv_diff_solve,'alpha1')
         alpha1= inv_model.ab_tv_diff_solve.alpha1;
      end
      if isfield(inv_model.ab_tv_diff_solve,'beta')
         beta= inv_model.ab_tv_diff_solve.beta;
      end
   end
   alpha2= calc_hyperparameter( inv_model);

   maxiter= 10;
   tol= 1e-4;
   keepiters= 0;
   if isfield(inv_model,'parameters')
      tol =     inv_model.parameters.term_tolerance;
      maxiter = inv_model.parameters.max_iterations;

      if isfield(inv_model.parameters,'keep_iterations')
         keepiters = inv_model.parameters.keep_iterations;
      end
   end
