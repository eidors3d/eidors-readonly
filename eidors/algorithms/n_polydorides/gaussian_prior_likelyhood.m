function likelihood= gaussian_prior_likelyhood( inv_model, x, y, J )
% Parameters for image
%   inv_model.gaussian_prior_likelihood.img_mean -> image mean
%   inv_model.gaussian_prior_likelihood.R_prior -> L*L' = inv(image covariance)
%   inv_model.gaussian_prior_likelihood.img_exp -> ( default = 2)
% Parameters for data
%   inv_model.gaussian_prior_likelihood.Noise -> L*L' = inv(Noise covariance)
%   inv_model.gaussian_prior_likelihood.data_exp -> ( default = 2)


x_m  = inv_model.gaussian_prior_likelihood.img_mean;
L_x  = inv_model.gaussian_prior_likelihood.R_prior;
p_x  = inv_model.gaussian_prior_likelihood.img_exp;
L_n  = inv_model.gaussian_prior_likelihood.Noise;
p_n  = inv_model.gaussian_prior_likelihood.data_exp;

img_residual=  x - x_m;
data_residual= y - J*x;

likelihood= exp(- norm(L_n * data_residual, p_n) ...
                - norm(L_x * img_residual,  p_x));