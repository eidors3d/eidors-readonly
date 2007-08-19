function jacobian = dg_calc_jacobian(img)
%% Solve the Jacobian matrix for a forward problem with a mapping function
%
% This function first applies a mapping function which gives the
% conductivities in each fine elements of the FEM structure from a set of
% parameters and auxiliary data:
%
% img.params_mapping.function = handle pointing on the mapping function
% img.params_mapping.params   = the parameters from which the conductivites
%                               are obtained
% img.params_mapping.perturb  = the perturbation applied to each parameters
%                               to compute the Jacobian
% imp.params_mapping.data = auxiliary data needed by the mapping function
%
% Dominique Gibert, April 2007
%
%%
[data_ref,img] = dg_fwd_solve(img);
n_data= size(data_ref.meas,1);
n_params= size(img.params_mapping.params,1);
jacobian= zeros(n_data,n_params);
for k= 1:n_params
    img.params_mapping.params(k)=img.params_mapping.params(k)+img.params_mapping.perturb(k);
    [data,img] = dg_fwd_solve(img);
    img.params_mapping.params(k)=img.params_mapping.params(k)-img.params_mapping.perturb(k);
    jacobian(:,k)= (data.meas-data_ref.meas);
end
