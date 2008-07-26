% Lung images
% $Id$

load montreal_data_1995
vh= zc_h_stomach_pre; % abdomen before fluid
vi= zc_stomach_0_5_60min(:,1); % right after drink

% GN solution - Gaussian prior
imdl.RtR_prior= @gaussian_HPF_prior;
imdl.solve=     @aa_inv_solve;
imdl.hyperparameter.value=1e-2;
img= inv_solve(imdl, vh, vi);

subplot(221); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl.RtR_prior= @noser_image_prior;
imdl.solve=     @aa_inv_solve;
imdl.hyperparameter.value=8e-2;
img= inv_solve(imdl, vh, vi);

subplot(222); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl= rmfield(imdl,'RtR_prior');
imdl.R_prior=   @ab_calc_tv_prior;
imdl.solve=     @ab_tv_diff_solve;
imdl.hyperparameter.value=1e-2;
img= inv_solve(imdl, vh, vi);

subplot(223); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl.R_prior=   @ab_calc_tv_prior;
imdl.solve=     @ab_tv_diff_solve;
imdl.hyperparameter.value=1e-1;
img= inv_solve(imdl, vh, vi);

subplot(224); show_fem(img); axis equal; axis off


print -r100 -dpng tutorial410d.png;
