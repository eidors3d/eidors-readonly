function img= ab_tv_diff_solve( inv_model, data1, data2)
% AB_TV_DIFF_SOLVE inverse solver for Andrea Borsic's
%   Total Variation solver for use with difference EIT
% img= ab_tv_diff_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: ab_tv_diff_solve.m,v 1.3 2005-12-02 15:28:46 aadler Exp $


hp= calc_hyperparameter( inv_model);
alpha1= hp(1); % Tikhonov parameter, should be large to start of right
alpha2= hp(2);
if isfield(inv_model,'parameters')
    maxiter = inv_model.parameters.max_iterations;
else
    maxiter= 3;
end
dva= data1 - data2;
%dva= data1;

sol=primaldual_tvrecon_lsearch(inv_model, dva ,maxiter,alpha1,alpha2);

img.name= 'solved by ab_tv_diff_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
