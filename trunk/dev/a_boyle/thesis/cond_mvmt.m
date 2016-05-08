% (C) 2016, Alistair Boyle
% License: GPL version 2 or version 3
imdlAMC = imdlA;
if isfield(imdlAMC, 'jacobian_bkgnd')
   imdlAMC = rmfield(imdlAMC, 'jacobian_bkgnd' );% allow jacobian_bkgnd estimation
end
if isfield(imdlAMC.fwd_model, 'jacobian')
   jacobian = imdlAMC.fwd_model.jacobian;
   imdlAMC.fwd_model = rmfield(imdlAMC.fwd_model, 'jacobian'); % replace jacobian with our own
else
   jacobian = imdlAMC.inv_solve_gn.jacobian{1};
end

imdlAMC.reconst_type= 'difference';
imdlAMC.solve = @inv_solve_gn;
nt = size(imdlAMC.fwd_model.elems,1);
ne = length(imdlAMC.fwd_model.electrode);
imdlAMC.hyperparameter.value = 1;
imdlAMC.inv_solve_gn.meas_input   = 'voltage';
imdlAMC.inv_solve_gn.meas_working = 'apparent_resistivity';
imdlAMC.inv_solve_gn.elem_prior   = {     'resistivity', 'movement'};
imdlAMC.inv_solve_gn.RtR_prior    = {  @prior_laplace  , @prior_movement_only_tikhonov };
imdlAMC.inv_solve_gn.elem_len     = {        nt        ,  ne       };
imdlAMC.inv_solve_gn.elem_fixed   = {        []        ,  [1:3 ne-2:ne]   }; % 3 first and 3 last electrodes are fixed
imdlAMC.inv_solve_gn.elem_working = {'log_conductivity', 'movement'};
imdlAMC.inv_solve_gn.elem_output  = {     'resistivity', 'movement'};
imdlAMC.inv_solve_gn.jacobian     = { jacobian, @jacobian_movement_only };
imdlAMC.inv_solve_gn.update_img_func = @mv_elec;
imdlAMC.inv_solve_gn.max_iterations = 10;
imdlAMC.inv_solve_gn.dtol           = -0.01; % stop at 1% improvement over the previous iteration
imdlAMC.inv_solve_gn.hyperparameter = {   1e-1         , [7e-2 5e-2*[1 1 1] 3e-2] };
imdlAMC.inv_solve_gn.return_working_variables = 1;
imdlAMC.inv_solve_gn.verbose=5; % > 8 = SVD... slow
imdlAMC.inv_solve_gn.meas_working = 'voltage';
imdlAMC.inv_solve_gn.fig_prefix = ['results/' name '-gn_mvmt'];
imdlAMC.inv_solve_gn.prior_data   = {    1  ,      0   };
imdlAMC.inv_solve_gn.prior_data{1} = imgAa.elem_data;
prior_move = (1:32)*0; fixed_elec = [1:3 30:32]; prior_move(fixed_elec) = duvw(fixed_elec,2)-0.04*(sel_line==1);
imdlAMC.inv_solve_gn.prior_data{2} = prior_move;

%%%% solve for conductivity and movement %%%%
imgAMCb = inv_solve( imdlAMC, va, vb);

   % THESIS figure!
dx = imgAMCb.elem_data(imgAMCb.params_sel{2}); dx = dx + prior_move(:);
clf; img=imgAMCb;
     img.elem_data = img.elem_data(img.params_sel{1});
     img.current_params = img.current_params{1}; img= rmfield(img, 'params_sel');
     img.calc_colours.clim = 50;
     img.calc_colours.ref_level = 0;
     show_fem_hollinhill(img);
   print('-dpdf', ['results/' name '-imgAMCbC.pdf']);
   eval(['! pdfcrop results/' name '-imgAMCbC.pdf']);
   eval(['! mv results/' name '-imgAMCbC-crop.pdf results/' name '-imgAMCbC.pdf']);
   eval(['! convert -verbose -density 200 -trim results/' name '-imgAMCbC.pdf -quality 90 results/' name '-imgAMCbC.jpg']);
clf; subplot(211); h=bar([ duvw(:,2)-0.04*(sel_line==1) (sum(dx,2))], 'EdgeColor', 'none');
        %hold on; plot([1:32]'*[1 1],[0.2*ones(1,32)]'*[1 -1] + duvw(:,2)*[1 1],'--','Color',[0.5 0.5 0.5]); hold off;
        hold on; plot([1:32]'*[1 1],[0.2*ones(1,32)]'*[1 -1],'--','Color',[0.5 0.5 0.5]); hold off;
        legend({'true movement','reconstructed'},'Location','SouthEast','Orientation','Horizontal','Box','off'); box off;
        ylabel('electrode mvmt [m]'); xlabel('electrode #'); axis tight;
        set(h(1),'FaceColor',[0.85 0.33 0.10]);
        set(h(2),'FaceColor',[0.47 0.67 0.19]);
        text(20,+0.35,'upslope'); text(20,-0.35,'downslope'); text(25,0.35,'0.2 m');
        ylim(common_ylim);
   print('-dpdf', ['results/' name '-imgAMCbM.pdf']);
   eval(['! pdfcrop results/' name '-imgAMCbM.pdf']);
   eval(['! mv results/' name '-imgAMCbM-crop.pdf results/' name '-imgAMCbM.pdf']);
