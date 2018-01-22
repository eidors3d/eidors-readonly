inv_crime=0;
%
% if(inv_crime==1)
%     % Make an ivnerse model - standard foward model inside
%     imdl = mk_common_model('c2C',16);
%     fmdl = imdl.fwd_model; %Extract model
%     fmdl = fix_model(fmdl);
%     imdl.fwd_model = fmdl;
%
%     %Stimulation - 16 elecs adjacent current/adjacent voltage
%     stim = mk_stim_patterns(16,1,'{ad}','{ad}');
%     fmdl.stimulation = stim; %Add to model
%     fmdl.approx_type='tri3';
%
%     cond_bkg = 1;
%     cond_obj = 2;
%
%
%     sim_img = mk_image(fmdl,cond_bkg);
%
%     %figure; show_fem(img,[0,1,3])
%     sim_img.fwd_solve.get_all_meas=1;
%     sim_img.fwd_model.stimulation=stim;
%
% %     sim_img.fwd_model.M_tensor.a = ones(size(sim_img.elem_data,1),1);
% %     sim_img.fwd_model.M_tensor.b = ones(size(sim_img.elem_data,1),1);
% %     sim_img.fwd_model.M_tensor.rot = zeros(size(sim_img.elem_data,1),1);
%
% %
%
%     homog_img=sim_img;
%
%     data_hom = fwd_solve(homog_img);
%     %pixel_group = [1,2,3,4];
%     %pixel_group = [115,138,95,137]; %b2C
%     pixel_group = [327,364,328,292,259,291]; %c2C
%
%     sim_img.elem_data(pixel_group) = cond_obj;
%
%     data = fwd_solve(sim_img)
%     data = add_noise(100, data, data_hom)
% else
cond_bkg = 1;
cond_obj = 2;


%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');

%Make a uniform image with unit background

fmdl= ng_mk_cyl_models(0,[16],[0.2,0,0.1]);
fmdl=fix_model(fmdl);
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

%Inverse model
imdl.solve='eidors_default';
%   imdl.hyperparameter=0.03;
imdl.RtR_prior='eidors_default';
imdl.jacobian_bkgnd.value=1;
imdl.reconst_type='difference';
imdl.fwd_model = fmdl;
imdl.name='built model';
imdl.type='inv_model';

%Notation of homog_img before
homog_img = mk_image(fmdl,cond_bkg);
homog_img.fwd_solve.get_all_meas=1;
homog_img.fwd_model.stimulation=stim;
homog_img.fwd_model.M_tensor.a = ones(size(homog_img.elem_data,1),1);
homog_img.fwd_model.M_tensor.b = ones(size(homog_img.elem_data,1),1);
homog_img.fwd_model.M_tensor.rot = zeros(size(homog_img.elem_data,1),1);

data_hom = fwd_solve(homog_img);


%% Ptensor options
fmdl = calc_closest_ellipse(fmdl);
% Default options
opts = [];
opts.H0_type = 'ptensor';
opts.use_hyper = 1;
opts.neumann = 'freespace';
opts.ptensor_its = 100;
opts.max_its = 100;
opts.update_delta = 1;
opts.inv_crime=inv_crime;
opts.update_U0 = 0;
opts.flexible = true;

opts.neumann = 'freespace';
opts.flexible = true;
opts.update_U0 = 1;


%% test scenarios
sepv=[0.3,0.4,0.5];
radv=[0.1,0.2,0.25];

nsep = length(sepv);
nrad = length(radv);

img_pt = cell(nsep, nrad);
r_pt = img_pt;
cts_pt = img_pt;
er_f11 = img_pt;
img_eid_diff = img_pt;
img_eid_abs = img_pt;

t_pt = zeros(nsep, nrad);
t_diff = t_pt;
t_eabs = t_pt;


for ii = 1:nsep
    for jj = 1:nrad
        
        % Sim data
        sep = sepv(ii);
        rad = radv(jj);
        sb1 = sprintf('solid ball1 = cylinder(%0.1f,%0.1f,0;%0.1f,%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;',sep,sep,sep,sep,rad);
        sb2 = sprintf('solid ball2 = cylinder(-%0.1f,-%0.1f,0;-%0.1f,-%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;',sep,sep,sep,sep,rad);
        extra={'ball1','ball2',[sb1,sb2]};
        
        fmdl_t= ng_mk_cyl_models(0.,[16],[0.2,0,0.05],extra);
        fmdl_t=fix_model(fmdl_t);
        fmdl_t.stimulation = stim; %Add to model
        fmdl_t.approx_type='tri3';
        sim_img= mk_image(fmdl_t, cond_bkg );
        sim_img.elem_data(fmdl_t.mat_idx{2})=cond_obj;
        %   img.elem_data(fmdl.mat_idx{3})=3;
        %   figure; show_fem(img);
        pixel_group=[fmdl_t.mat_idx{2}];
        %pixel_group = [102,123,103,83,66,82]; %b2C
        %pixel_group = [327,364,328,292,259,291]; %c2C
        %pixel_group = 1:length(fmdl.elems(:,1));
        
        data = fwd_solve(sim_img);
        data = add_noise(100, data, data_hom);
        
        
        %% Eidors in-built inverse diff solver
        imdl.hyperparameter.value = 0.0003;
        imdl.reconst_type = 'difference';
        tic
        img_eid_diff{ii,jj} = inv_solve(imdl, data_hom,data);
        t_diff = toc;
%         
%         figure(1);
%         subplot(121); show_fem(sim_img,1)
%         subplot(122); show_fem(img_eid_diff,1)
        
        %% Eidors in-built inverse abs solver
        imdl.fwd_model=fmdl;
        imdl.solve =  @inv_solve_core;
        imdl.reconst_type = 'absolute';
        
        tic
        img_eid_abs{ii,jj} = inv_solve(imdl, data);
        t_eabs = toc;
        
%         figure(2);
%         subplot(121); show_fem(sim_img,1)
%         subplot(122); show_fem(img_eid_abs,1)
%         
%         data_rec=fwd_solve(img_eid_abs);
%         figure(3); plot(data.meas,'r'); hold on; plot(data_hom.meas,'b');  hold on; plot(data_rec.meas,'b')

        
        %% P-Tensor lbfgs solvers
        
        tic
        [img_pt{ii,jj}, r_pt{ii,jj}, cts_pt{ii,jj}, er_f11{ii,jj}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
        t_pt(ii,jj) = toc;
        hold all
        
        
    end
end













% figure(4)
% 
% % % 3 update settings in freespace
% % opts.neumann = 'freespace';
% % opts.update_U0 = 0;
% % opts.flexible = false;
% % tic
% % [x_f00, r_f00, cts_f00, er_f00] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% % tf00 = toc
% % semilogy(0:100, r_f00/ r_f00(1))
% 
% %
% 
% opts.neumann = 'freespace';
% 
% 
% hold all
% semilogy(0:100, r_f10/ r_f10(1))
% 
% %
% opts.neumann = 'freespace';
% opts.flexible = true;
% opts.update_U0 = 1;
% tic
% [x_f11, r_f11, cts_f11, er_f11] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% tf11 = toc
% hold all
% semilogy(0:100, r_f11/ r_f11(1))
% 
% % %
% % opts.ptensor_its = 20;
% % [x_f20, r_f20, cts_f20] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% % semilogy(0:50, r_f20/ r_f0(1))
% 
% %%
% % % 3 update settings on a disc
% % opts.neumann = 'disc';
% % opts.update_U0 = 0;
% % opts.flexible = false;
% % tic
% % [x_d00, r_d00, cts_d00, er_d00] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% % td00 = toc
% % semilogy(0:100, r_d00/ r_f00(1))
% %
% % %
% % opts.neumann = 'disc';
% % opts.update_U0 = 0;
% % opts.flexible = true;
% % tic
% % [x_d10, r_d10, cts_d10, er_d10] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% % td10 = toc
% % semilogy(0:100, r_d10/ r_f00(1))
% %
% % %
% % opts.neumann = 'disc';
% % opts.update_U0 = 1;
% % opts.flexible = true;
% % tic
% % [x_d11, r_d11, cts_d11, er_d11] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% % td11 = toc
% % semilogy(0:100, r_d11/ r_f00(1))
% 
% %%
% opts.neumann = 'freespace';
% opts.update_U0 = 1;
% opts.flexible = true;
% opts.H0_type = 'ptensor_full';
% tic;
% [x_pf, r_pf, cts_pf, er_pf] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% t_pf = toc
% semilogy(0:100, r_pf/r_pf(1))
% % %
% % opts.ptensor_its = 20;
% % [x_d20, r_d20] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% % semilogy(0:100, r_d20/ r_f0(1))
% 
% 
% %
% % compare to DGN0
% 
% 
% 
% % % Include reg
% % opts.update_U0 = 0;
% % opts.use_hyper = 1;
% % [x_Phy, r_Phy, cts_Phy] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% %
% %
% %
% % % Include reg and recalc U0
% % opts.update_U0 = 1;
% % opts.use_hyper = 1;
% % [x_Phyup, r_Phyup, cts_Phyup] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% %
% % % Re-scaling
% % opts.use_hyper = 0;
% % opts.rescale = 1;
% % [x_Pre, r_Pre, cts_Pre] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% 
% 
% %%
% % % no ptensor
% % opts.H0_type = '';
% % tic
% % [x_I, r_I, cts_I, er_I] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% % tI = toc
% % semilogy(0:100, r_I/ r_f00(1))
% 
% %%
% opts.H0_type = 'DGN0';
% opts.flexible = false;
% tic
% [x_GN0, r_GN0, cts_GN0, er_GN0] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% tGN0 = toc
% semilogy(0:100, r_GN0/ r_GN0(1))
% 
% %%
% opts.H0_type = 'DGN0';
% opts.flexible = true;
% tic
% [x_GN1, r_GN1, cts_GN1, er_GN1] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% tGN1 = toc
% semilogy(0:100, r_GN1/ r_GN1(1))
% 
% %%
% opts.H0_type = 'true';
% opts.flexible = false;
% tic
% [x_t0, r_t0, cts_t0, er_t0] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% tt0 = toc
% semilogy(0:100, r_t0/ r_t0(1))
% 
% %
% opts.H0_type = 'true';
% opts.flexible = true;
% tic
% [x_t1, r_t1, cts_t1, er_t1] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img);
% tt1 = toc
% semilogy(0:100, r_t1/ r_t1(1))
% 
% %%
% hold off
% figure(5)
% % semilogy(0:100, er_f00/ er_f00(1))
% hold all
% semilogy(0:100, er_f10/ er_f10(1))
% semilogy(0:100, er_f11/ er_f11(1))
% % semilogy(0:100, er_d00/ er_d00(1))
% % semilogy(0:100, er_d10/ er_d10(1))
% % semilogy(0:100, er_d11/ er_d11(1))
% semilogy(0:100, er_pf/er_pf(1))
% semilogy(0:100, er_GN0/ er_GN0(1))
% semilogy(0:100, er_GN1/ er_GN1(1))
% semilogy(0:100, er_t0/ er_t0(1))
% semilogy(0:100, er_t1/ er_t1(1))
