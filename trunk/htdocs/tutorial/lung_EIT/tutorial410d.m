% Abdomen Images $Id$

load montreal_data_1995
vh= zc_h_stomach_pre; % abdomen before fluid
vi= zc_stomach_0_5_60min(:,1); % right after drink

% GN solution - Gaussian prior
imdl.RtR_prior= @prior_gaussian_HPF;
imdl.solve=     @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value=1e-2;
img= inv_solve(imdl, vh, vi);

subplot(221); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl.RtR_prior= @prior_noser;
imdl.solve=     @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value=8e-2;
img= inv_solve(imdl, vh, vi);

subplot(222); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl= rmfield(imdl,'RtR_prior');
imdl.R_prior=   @prior_TV;
imdl.solve=     @inv_solve_TV_pdipm;
imdl.hyperparameter.value=3e-4;
img= inv_solve(imdl, vh, vi);

subplot(223); show_fem(img); axis equal; axis off

% GN solution - Noser prior
imdl.R_prior=   @prior_TV;
imdl.solve=     @inv_solve_TV_pdipm;
imdl.hyperparameter.value=3e-3;
img= inv_solve(imdl, vh, vi);

subplot(224); show_fem(img); axis equal; axis off


print_convert tutorial410d.png '-density 150';
