function img= ab_tv_diff_solve( inv_model, data1, data2)
% AB_TV_DIFF_SOLVE inverse solver for Andrea Borsic's
%   Total Variation solver for use with difference EIT
% img= ab_tv_diff_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
% Parameters
%   alpha1
%   alpha2

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$


[alpha1,alpha2,beta,maxiter,tol,keepiters]= get_params(inv_model);

dva = calc_difference_data( data1, data2, inv_model.fwd_model);
% TEST CODE -> Put elsewhere
back_val = get_good_background(inv_model, data1);
inv_model.jacobian_bkgnd.value= back_val;

sol= [];
for i=1:size(dva,2)
   soln=primaldual_tvrecon_lsearch(inv_model, dva(:,i), ...
        maxiter,alpha1,alpha2, tol, beta);

   if ~keepiters
      soln=soln(:,end);
   end

   sol=[sol, soln];
end

img.name= 'solved by ab_tv_diff_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;

function [alpha1,alpha2,beta,maxiter,tol,keepiters]= ...
          get_params(inv_model);
   alpha1= 1e-2;
   beta= 1e-4;
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
   try 
      tol =     inv_model.parameters.term_tolerance;
   end
   try
      maxiter = inv_model.parameters.max_iterations;
   end
   try
      keepiters = inv_model.parameters.keep_iterations;
   end

function back_val = get_good_background(inv_mdl, data1);

   % Create homogeneous model
   IM= eidors_obj('image','');
   IM.fwd_model= inv_mdl.fwd_model;
   s= ones(size(IM.fwd_model.elems,1),1);
   IM.elem_data= s;

   vsim= fwd_solve( IM);
   back_val=abs( data1\vsim.meas );
   back_val=1;

