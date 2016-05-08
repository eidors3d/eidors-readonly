% (C) 2016, Alistair Boyle
% License: GPL version 2 or version 3
if isfield(imdlA, 'jacobian_bkgnd')
   imdlA = rmfield(imdlA, 'jacobian_bkgnd' );% allow jacobian_bkgnd estimation
end
if isfield(imdlA.fwd_model, 'jacobian')
   jacobian = imdlA.fwd_model.jacobian;
   imdlA.fwd_model = rmfield(imdlA.fwd_model, 'jacobian'); % replace jacobian with our own
else
   jacobian = imdlA.inv_solve_gn.jacobian{1};
end
imdlA.reconst_type= 'absolute';
imdlA.solve = @inv_solve_gn;
nt = size(imdlA.fwd_model.elems,1);
ne = length(imdlA.fwd_model.electrode);
imdlA.hyperparameter.value = 1;
imdlA.inv_solve_gn.meas_input   = 'voltage';
imdlA.inv_solve_gn.meas_working = 'apparent_resistivity';
imdlA.inv_solve_gn.elem_prior   = 'resistivity';
imdlA.inv_solve_gn.RtR_prior    = @prior_laplace;
imdlA.inv_solve_gn.elem_len     = nt;
imdlA.inv_solve_gn.elem_fixed   = [];
imdlA.inv_solve_gn.elem_working = 'log_conductivity';
imdlA.inv_solve_gn.elem_output  = 'log10_resistivity';
imdlA.inv_solve_gn.jacobian     = jacobian;
imdlA.inv_solve_gn.max_iterations = 10;
imdlA.inv_solve_gn.dtol           = -0.01; % stop at 1% improvement over the previous iteration
imdlA.inv_solve_gn.return_working_variables = 1;
imdlA.inv_solve_gn.verbose=5; % > 8 = SVD... slow
imdlA.inv_solve_gn.show_fem = @ab_show_fem;
imdlA.inv_solve_gn.hyperparameter = [1.25e2 1.25e1 1.25e0];
imdlA.inv_solve_gn.fig_prefix = ['results/' name '-Aa'];
imdlA.inv_solve_gn.prior_data = BACKGROUND_RAa;

img= inv_solve( imdlA, va);

% THESIS figure
imgAa.calc_colours.clim = 50;
imgAa.calc_colours.ref_level = 50;
clf; show_fem_hollinhill(imgAa); drawnow;
print('-dpdf', ['results/' name '-imgAa_log10.pdf']);
eval(['! pdfcrop results/' name '-imgAa_log10.pdf']);
eval(['! mv      results/' name '-imgAa_log10-crop.pdf results/' name '-imgAa_log10.pdf']);
eval(['! convert -verbose -density 200 -trim results/' name '-imgAa_log10.pdf -quality 90 results/' name '-imgAa_log10.jpg']);
% THESIS figure
clf; show_fem_hollinhill(imgAa,0); drawnow;
print('-dpdf', ['results/' name '-imgAa_log10_full.pdf']);
eval(['! pdfcrop results/' name '-imgAa_log10_full.pdf']);
eval(['! mv      results/' name '-imgAa_log10_full-crop.pdf results/' name '-imgAa_log10_full.pdf']);
eval(['! convert -verbose -density 200 -trim results/' name '-imgAa_log10_full.pdf -quality 90 results/' name '-imgAa_log10_full.jpg']);
