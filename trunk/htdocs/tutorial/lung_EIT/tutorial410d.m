% Lung images
% $Id$

load montreal_data_1995
vh= zc_h_stomach_pre; % abdomen before fluid
vi= zc_stomach_0_5_60min(:,1); % right after drink

% GN solution - Gaussian prior
imdl.RtR_prior= @gaussian_HPF_prior;
imdl.solve=     @GN_one_step_diff_solve;
imdl.hyperparameter.value=1e-2;
img= inv_solve(imdl, vh, vi);

subplot(221); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl.RtR_prior= @noser_image_prior;
imdl.solve=     @GN_one_step_diff_solve;
imdl.hyperparameter.value=8e-2;
img= inv_solve(imdl, vh, vi);

subplot(222); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl= rmfield(imdl,'RtR_prior');
imdl.R_prior=   @calc_TV_prior;
imdl.solve=     @TV_diffusivity_solve;
imdl.hyperparameter.value=3e-4;
img= inv_solve(imdl, vh, vi);

subplot(223); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl.R_prior=   @calc_TV_prior;
imdl.solve=     @TV_diffusivity_solve;
imdl.hyperparameter.value=3e-3;
img= inv_solve(imdl, vh, vi);

subplot(224); show_fem(img); axis equal; axis off


print_convert tutorial410d.png '-density 150';
