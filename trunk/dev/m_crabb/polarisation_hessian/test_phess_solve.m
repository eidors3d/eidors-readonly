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
%     sim_img{ii,jj} = mk_image(fmdl,cond_bkg);
%
%     %figure; show_fem(img,[0,1,3])
%     sim_img{ii,jj}.fwd_solve.get_all_meas=1;
%     sim_img{ii,jj}.fwd_model.stimulation=stim;
%
% %     sim_img{ii,jj}.fwd_model.M_tensor.a = ones(size(sim_img{ii,jj}.elem_data,1),1);
% %     sim_img{ii,jj}.fwd_model.M_tensor.b = ones(size(sim_img{ii,jj}.elem_data,1),1);
% %     sim_img{ii,jj}.fwd_model.M_tensor.rot = zeros(size(sim_img{ii,jj}.elem_data,1),1);
%
% %
%
%     homog_img=sim_img{ii,jj};
%
%     data_hom = fwd_solve(homog_img);
%     %pixel_group = [1,2,3,4];
%     %pixel_group = [115,138,95,137]; %b2C
%     pixel_group = [327,364,328,292,259,291]; %c2C
%
%     sim_img{ii,jj}.elem_data(pixel_group) = cond_obj;
%
%     data = fwd_solve(sim_img{ii,jj})
%     data = add_noise(100, data, data_hom)
% else



%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');


%Make a uniform image with unit background
cond_bkg = 1;

fmdl= ng_mk_cyl_models([0,1,0.1],[16],[0.2,0,0.1]);
%fmdl= ng_mk_cyl_models(0,[16],[0.2,0,0.1]);

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
% homog_img.fwd_model.M_tensor.a = ones(size(homog_img.elem_data,1),1);
% homog_img.fwd_model.M_tensor.b = ones(size(homog_img.elem_data,1),1);
% homog_img.fwd_model.M_tensor.rot = zeros(size(homog_img.elem_data,1),1);
fmdl = calc_closest_ellipse(fmdl);

data_hom = fwd_solve(homog_img);


%% Ptensor options
fmdl = calc_closest_ellipse(fmdl);
% Default options
% popts = [];
% popts.H0_type = 'ptensor';
% popts.use_hyper = 1;
% popts.neumann = 'disc';
% popts.ptensor_its = 100;
% popts.max_its = 100;
% popts.update_delta = 1;
% popts.inv_crime=inv_crime;
% popts.update_U0 = 1;
% popts.flexible = true;

%opts.neumann = 'disc';
popts.neumann = 'freespace';

popts.flexible = true;
popts.update_U0 = 1;


%% test scenarios
sepv=0.25;%logspace(log10(0.25), log10(0.35),2);
radv=0.25;%logspace(log10(0.16), log10(0.25),2);
conv=0.8;%linspace(0.8,2,3);
%%
% sepv=0.2
% radv =0.1646
% conv=1.5

cond_obj1 = 2;
cond_obj2 = conv;%1.5;

nsep = length(sepv);
nrad = length(radv);
ncon = length(conv);

img_pt = cell(nsep, nrad);
r_pt = img_pt;
er_f11 = img_pt;
img_eid_diff = img_pt;
img_eid_abs = img_pt;

t_pt = zeros(nsep, nrad);
t_diff = t_pt;
t_eabs = t_pt;
eres = t_pt;
cts_pt = t_pt;


for ii = 1:nsep
    for jj = 1:nrad
        for    kk=1:ncon
            
            
            % Sim data
            sep = sepv(ii);
            rad = radv(jj);
            sb1 = sprintf('solid ball1 = cylinder(%0.1f,%0.1f,0;%0.1f,%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;',sep,sep,sep,sep,rad);
            sb2 = sprintf('solid ball2 = cylinder(-%0.1f,-%0.1f,0;-%0.1f,-%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) and not cylinder(%0.1f,%0.1f,0;%0.1f,%0.1f,1;%0.1f) -maxh=0.1;',sep,sep,sep,sep,rad, sep,sep,sep,sep,rad);
            extra={'ball1','ball2',[sb1,sb2]};

            fmdl_t= ng_mk_cyl_models(0.,[16],[0.2,0,0.05],extra);
            fmdl_t=fix_model(fmdl_t);
            fmdl_t.stimulation = stim; %Add to model
            fmdl_t.approx_type='tri3';
            sim_img{ii,jj}= mk_image(fmdl_t, cond_bkg );
            sim_img{ii,jj}.elem_data(fmdl_t.mat_idx{2})=cond_obj1;
            sim_img{ii,jj}.elem_data(fmdl_t.mat_idx{3})=cond_obj2(kk);
            %   figure; show_fem(img);
            pixel_group=[fmdl_t.mat_idx{2}];
            %pixel_group = [102,123,103,83,66,82]; %b2C
            %pixel_group = [327,364,328,292,259,291]; %c2C
            %pixel_group = 1:length(fmdl.elems(:,1));

            n_level = 50;

            data = fwd_solve(sim_img{ii,jj});
%             figure; plot(data.meas,'b'); hold on;
            data = add_noise(n_level, data, data_hom);
%             plot(data.meas,'r');

            %% Eidors in-built inverse diff solver

            h_param_eid=0.00005;


            imdl.hyperparameter.value =   h_param_eid;
            imdl.reconst_type = 'difference';
            tic
            img_eid_diff{ii,jj} = inv_solve(imdl, data_hom,data);
            t_diff(ii,jj) = toc;
    %         
    %         figure(1);
    %         subplot(121); show_fem(sim_img{ii,jj},1)
    %         subplot(122); show_fem(img_eid_diff,1)

            %% Eidors in-built inverse abs solver
            imdl.fwd_model=fmdl;
            imdl.hyperparameter.value = h_param_eid ;

            imdl.solve =  @inv_solve_core;
            imdl.reconst_type = 'absolute';

            opt = [];
    %         opt.
            imdl.inv_solve_core.verbose=0;
            imdl.inv_solve_core.tol=10^-12;


            tic
            img_eid_abs{ii,jj} = inv_solve(imdl, data);
            t_eabs(ii,jj) = toc;

    %         figure(2);
    %         subplot(121); show_fem(sim_img{ii,jj},1)
    %         subplot(122); show_fem(img_eid_abs,1)
    %         
    %         data_rec=fwd_solve(img_eid_abs);
    %         figure(3); plot(data.meas,'r'); hold on; plot(data_hom.meas,'b');  hold on; plot(data_rec.meas,'b')

            %% P-Tensor lbfgs solvers        
            %p-tensor jacobian approx factor 2 out
            imdl.hyperparameter.value = 2*h_param_eid;
            popts = [];
            popts.d_tol = 1e-4;
            popts.c1 = 1e-2;
            popts.c2 = 0.75;
            popts.max_its = 100;
            popts.use_hyper = 0;
            popts.update_U0 = 1;
            popts.H0_type = 'ptensor';
            popts.neumann = 'freespace';
            
%             tic
%             [img_pt_free{ii,jj}, r_pt_free{ii,jj}, cts_pt_free(ii,jj), er_f11_free{ii,jj}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, popts, sim_img{ii,jj});
%             t_pt_free(ii,jj) = toc;
%             eres_free(ii,jj) = r_pt_free{ii,jj}(end);
            
            %% P-Tensor lbfgs solvers
            imdl.hyperparameter.value = 2*h_param_eid;
            popts.neumann = 'disc';
            
            tic
            [img_pt_disc{ii,jj}, r_pt_disc{ii,jj}, cts_pt_disc(ii,jj), er_f11_disc{ii,jj}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, popts, sim_img{ii,jj});
            t_pt_disc(ii,jj) = toc;
            eres_disc(ii,jj) = r_pt_disc{ii,jj}(end);
           
           
        
            %% GN lbfgs solver
            imdl.hyperparameter.value = 2*h_param_eid;
            popts.H0_type = 'DGN0';
            popts.use_hyper = 1;
            
            tic
            [img_gn0{ii,jj}, r_gn0{ii,jj}, cts_gn0(ii,jj), er_gn0{ii,jj}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, popts, sim_img{ii,jj});
            t_gn0(ii,jj) = toc;
            eres_gn0(ii,jj) = r_gn0{ii,jj}(end);
            
            %% lambda I + hp^2 RtR lbfgs solver
            popts.H0_type = 'identity';
            popts.c1 = 1e-4;
            popts.c2 = 0.9;
             tic
            [img_I{ii,jj}, r_I{ii,jj}, cts_I(ii,jj), er_I{ii,jj}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, popts, sim_img{ii,jj});
            t_I(ii,jj) = toc;
            eres_I(ii,jj) = r_I{ii,jj}(end);

            

        end
    end
end





%% Show
for ii = 1:nsep
    for jj = 1:nrad
        for    kk=1:ncon
            
            figure
            subplot(2,2,1); show_fem(sim_img{ii,jj},1)
%             subplot(2,3,2); show_fem(img_eid_diff{ii,jj},1)
            subplot(2,2,2); show_fem(img_eid_abs{ii,jj},1)
            
            subplot(2,2,3); show_fem(img_pt_disc{ii,jj},1)            
%             subplot(2,3,3); show_fem(img_pt_free{ii,jj},1)
            subplot(2,2,4); show_fem(img_gn0{ii,jj},1);
            


        end
    end
end

% Separate loops as subplots were going all over the place ?
for ii = 1:nsep
    for jj = 1:nrad
        for    kk=1:ncon
            
            figure
%             semilogy(0:length(r_pt_free{ii,jj})-1,r_pt_free{ii,jj}/r_pt_free{ii,jj}(1))
            semilogy(0:length(r_pt_disc{ii,jj})-1,r_pt_disc{ii,jj}/r_pt_disc{ii,jj}(1))
            hold all

            semilogy(0:length(r_gn0{ii,jj})-1,r_gn0{ii,jj}/r_gn0{ii,jj}(1))
            semilogy(0:length(r_I{ii,jj})-1,r_I{ii,jj}/r_I{ii,jj}(1))
            grid on
            hold off


        end
    end
end


% figure(4)
% 
% % % 3 update settings in freespace
% % opts.neumann = 'freespace';
% % opts.update_U0 = 0;
% % opts.flexible = false;
% % tic
% % [x_f00, r_f00, cts_f00, er_f00] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
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
% [x_f11, r_f11, cts_f11, er_f11] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
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
% % [x_d00, r_d00, cts_d00, er_d00] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% % td00 = toc
% % semilogy(0:100, r_d00/ r_f00(1))
% %
% % %
% % opts.neumann = 'disc';
% % opts.update_U0 = 0;
% % opts.flexible = true;
% % tic
% % [x_d10, r_d10, cts_d10, er_d10] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% % td10 = toc
% % semilogy(0:100, r_d10/ r_f00(1))
% %
% % %
% % opts.neumann = 'disc';
% % opts.update_U0 = 1;
% % opts.flexible = true;
% % tic
% % [x_d11, r_d11, cts_d11, er_d11] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% % td11 = toc
% % semilogy(0:100, r_d11/ r_f00(1))
% 
% %%
% opts.neumann = 'freespace';
% opts.update_U0 = 1;
% opts.flexible = true;
% opts.H0_type = 'ptensor_full';
% tic;
% [x_pf, r_pf, cts_pf, er_pf] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
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
% % [x_I, r_I, cts_I, er_I] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% % tI = toc
% % semilogy(0:100, r_I/ r_f00(1))
% 
% %%
% opts.H0_type = 'DGN0';
% opts.flexible = false;
% tic
% [x_GN0, r_GN0, cts_GN0, er_GN0] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% tGN0 = toc
% semilogy(0:100, r_GN0/ r_GN0(1))
% 
% %%
% opts.H0_type = 'DGN0';
% opts.flexible = true;
% tic
% [x_GN1, r_GN1, cts_GN1, er_GN1] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% tGN1 = toc
% semilogy(0:100, r_GN1/ r_GN1(1))
% 
% %%
% opts.H0_type = 'true';
% opts.flexible = false;
% tic
% [x_t0, r_t0, cts_t0, er_t0] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
% tt0 = toc
% semilogy(0:100, r_t0/ r_t0(1))
% 
% %
% opts.H0_type = 'true';
% opts.flexible = true;
% tic
% [x_t1, r_t1, cts_t1, er_t1] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts, sim_img{ii,jj});
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
