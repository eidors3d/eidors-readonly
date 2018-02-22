% Compare eigenvectors/singular vectors of Jacobian and Hess approx to true
inv_crime = 0;


%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');


%Make a uniform image with unit background
cond_bkg = 1;

fmdl= ng_mk_cyl_models([0,1,0.1],[16],[0.2,0,0.1]);
% fmdl= ng_mk_cyl_models(0,[16],[0.2,0,0.1]);
fmdl=fix_model(fmdl);
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

%Inverse model
imdl.solve='eidors_default';
imdl.hyperparameter.value=1e-4;
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
opts = [];
opts.H0_type = 'ptensor';
opts.use_hyper = 1;
opts.neumann = 'freespace';
opts.ptensor_its = 100;
opts.max_its = 100;
opts.update_delta = 1;
opts.inv_crime=inv_crime;
opts.update_U0 = 1;
opts.flexible = true;

opts.neumann = 'freespace';
opts.flexible = true;
opts.update_U0 = 1;


%% test scenarios
sepv=0.3;%logspace(log10(0.22), log10(0.3),3);
radv=0.25;%logspace(log10(0.2), log10(0.35),3);


cond_obj1 = 2;
cond_obj2 = 1.5;

nsep = length(sepv);
nrad = length(radv);

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


ii=1; jj=1;

% for ii = 1:nsep
%     for jj = 1:nrad
        
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
        sim_img= mk_image(fmdl_t, cond_bkg );
        sim_img.elem_data(fmdl_t.mat_idx{2})=cond_obj1;
        sim_img.elem_data(fmdl_t.mat_idx{3})=cond_obj2;
        %   figure; show_fem(img);
%         pixel_group=[fmdl_t.mat_idx{2}];
        %pixel_group = [102,123,103,83,66,82]; %b2C
        %pixel_group = [327,364,328,292,259,291]; %c2C
        %pixel_group = 1:length(fmdl.elems(:,1));
        
        data = fwd_solve(sim_img);
        data = add_noise(100, data, data_hom);
        delta_d = data.meas - data_hom.meas;
        
        
        % True Jacobian & Hess fwd_model,img,elem_list,delta_d
        pixel_group = 1:size(fmdl.elems,1);
        J=calc_jacobian(homog_img);
        homog_img.fwd_model.approx_type='tri3';
        H = calc_hessian_obj(homog_img.fwd_model,homog_img,pixel_group, delta_d);
        H_diag = calc_hessian_diag(homog_img.fwd_model,homog_img,pixel_group);
        H = H + calc_hyperparameter(imdl)^2*calc_RtR_prior(imdl);
        
        %%
        
        % P-tensor approximations from free-space
        homog_img.fwd_solve.get_all_meas=1;
        u0= data_hom.volt;
        DU0 = calc_grad_potential(homog_img, u0);
        [Hii, du2, d2u, J_phess] = calc_phessian_obj(homog_img.fwd_model,homog_img,DU0,data_hom.meas,'disc');
        
        % Compare principal angles between subspaces
        % Choose threshold of ratio of singular values not in nullspace
        indx = 1:20;
        
        % Compare some singular vectors
        [Uj, Sj, Vj] = svd(J);
        [Upj, Spj, Vpj] = svd(J_phess);
        
        
        
        iter_solve = [5, 10, 15, 25, 50];
        for kk=1:length(iter_solve)
            
            % P-Tensor approx after X iterations of BFGS
            opt.d_tol = 0;
            opt.max_its = iter_solve(kk);
            opt.mem = iter_solve(kk);
            opt.use_hyper = 1;
            opt.H0_type = 'ptensor';
            [~,~,~,~,H_PBFGS{kk}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
            opt.H0_type = 'identity';
            [~,~,~,~,H_IBFGS{kk}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
            opt.H0_type = 'regu';
            [~,~,~,~,H_RBFGS{kk}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
            opt.H0_type = 'DGN0';
            [~,~,~,~,H_GNBFGS{kk}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
            
            
            
            
            % Reg contn
            RtR = imdl.hyperparameter.value^2*calc_RtR_prior(imdl);
            
            
            
            %         %
            %         H_PBFGS = H_vers{1};% - RtR;
            %         H_IBFGS = H_vers{2};
            %         H_RBFGS = H_vers{3} - RtR;
            
            [~, Shp{kk}, Vhp] = svd(H_PBFGS{kk});
            [~, Shi{kk}, Vhi] = svd(H_IBFGS{kk});
            [~, Shr{kk}, Vhr] = svd(H_RBFGS{kk});
            [~, Shg{kk}, Vhg] = svd(H_GNBFGS{kk});
            [~, Sh{kk}, Vh] = svd(H);
            
           
            
           
            
            
            
            % Hess prin angles
%             indx = diag(Sh{kk})/Sh{kk}(1,1) > thresh;
            H_true_basis = Vh(:,indx);
            
%             indx = diag(Shp{kk})/Shp{kk}(1,1) > thresh;
            H_p_basis = Vhp(:,indx);
            
%             indx = diag(Shi{kk})/Shi{kk}(1,1) > thresh;
            H_i_basis = Vhi(:,indx);
            
%             indx = diag(Shr{kk})/Shr{kk}(1,1) > thresh;
            H_r_basis = Vhr(:,indx);
            
            H_g_basis = Vhg(:,indx);
            
            % Prin angles by SVD
            cos_theta_Hp(:,kk) = svd(H_true_basis.'*H_p_basis);
            cos_theta_Hr(:,kk) = svd(H_true_basis.'*H_r_basis);
            cos_theta_Hi(:,kk) = svd(H_true_basis.'*H_i_basis);
            cos_theta_Hg(:,kk) = svd(H_true_basis.'*H_g_basis);
            
        end
        
        %%
        
        % Jacobian principal angles
        thresh = 0;
%         indx = diag(Spj)/Spj(1,1) > thresh;
        J_true_basis = Vj(:, indx);
        
%         indx = diag(Spj)/Spj(1,1) > thresh;
        J_p_basis = Vpj(:,indx);
        
        % Prin angles by SVD
        cos_theta_J = svd(J_true_basis.' * J_p_basis);

        
  
        
        
        
%     end
% end


%% Plot some things
figure(1)
plot(acos(cos_theta_Hp))
figure(2)
plot(acos(cos_theta_Hr))
figure(3)
plot(acos(cos_theta_Hi))

% figure(4)
% plot(acos([cos_theta_Hp(:,end), cos_theta_Hr(:,end), cos_theta_Hi(:,end)]))
% 
figure(5)
plot(acos(cos_theta_J))

%% Convergence to true Hess
iter_solve = [1:25];
for kk=1:length(iter_solve)
    
    % P-Tensor approx after X iterations of BFGS
    opt.d_tol = 0;
    opt.max_its = iter_solve(kk);
    opt.mem = iter_solve(kk);
    opt.use_hyper = 1;
    opt.H0_type = 'ptensor';
    [~,~,~,~,H_PBFGS{kk}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
    opt.H0_type = 'DGN0';
    [~,~,~,~,H_GNBFGS{kk}] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
    
    dHP(kk) = norm(H_PBFGS{kk} - H)/norm(H);
    dHG(kk) = norm(H_GNBFGS{kk} - H)/norm(H);

    
    
end

