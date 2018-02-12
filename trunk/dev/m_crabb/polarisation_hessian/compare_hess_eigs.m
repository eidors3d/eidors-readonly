% Compare eigenvectors/singular vectors of Jacobian and Hess approx to true
inv_crime = 0;


%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');


%Make a uniform image with unit background
cond_bkg = 1;

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


for ii = 1:nsep
    for jj = 1:nrad
        
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
        
        % P-tensor approximations from free-space
        homog_img.fwd_solve.get_all_meas=1;
        u0= data_hom.volt;
        DU0 = calc_grad_potential(homog_img, u0);
        [Hii, du2, d2u, J_phess] = calc_phessian_obj(homog_img.fwd_model,homog_img,DU0,data_hom.meas,'disc');
        
        
        % P-Tensor approx after 20 iterations of BFGS
        opt.max_its = 20;
        opt.mem = 20;
        [~,~,~,~,H_PBFGS, H_IBFGS] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opt, sim_img);
        
        % Reg contn
        RtR = calc_RtR_prior(imdl);
        
        % Compare some singular vectors
        [Uj, Sj, Vj] = svd(J);
        [Upj, Spj, Vpj] = svd(J_phess);
        
        %
        

        
    end
end