% Calculate Hessian of objective, for a fixed dataset, as an element's
% contrast is varied.  See also comparess_obj_saturation, for varied
% dataset on fixed homogeneous domain.


% % Data simulation model
% smdl = mk_common_model('c2C',16);
% fmdl = fix_model(smdl.fwd_model);
% stim = mk_stim_patterns(16,1,'{ad}','{ad}');
% fmdl.stim = stim;
% fmdl.approx_type = 'tri3';
% 
% smdl.fwd_model = fmdl;
% 
% % Inversion model
% imdl = mk_common_model('b2C',16);
% fmdl = fix_model(imdl.fwd_model);
% stim = mk_stim_patterns(16,1,'{ad}','{ad}');
% fmdl.stim = stim;
% fmdl.approx_type = 'tri3';
% 
% 
% imdl.fwd_model = fmdl;
% 
% 
% % homog_imgs
% cond_bkg = 1;
% 
% % Homog coarse homog_img
% homog_img = mk_image(imdl, cond_bkg);
% homog_img.fwd_model.M_tensor.a = ones(size(homog_img.elem_data,1),1);
% homog_img.fwd_model.M_tensor.b = ones(size(homog_img.elem_data,1),1);
% homog_img.fwd_model.M_tensor.rot = zeros(size(homog_img.elem_data,1),1);
% 
% homog_img = homog_img;
% 
% % Fine homog_img
% sim_homog_img = mk_image(smdl, cond_bkg);
% sim_homog_img.fwd_model.M_tensor.a = ones(size(sim_homog_img.elem_data,1),1);
% sim_homog_img.fwd_model.M_tensor.b = ones(size(sim_homog_img.elem_data,1),1);
% sim_homog_img.fwd_model.M_tensor.rot = zeros(size(sim_homog_img.elem_data,1),1);
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


cond_obj = 0.8:0.1:1.8;


nsep = length(sepv);
nrad = length(radv);

homog_img_pt = cell(nsep, nrad);
r_pt = homog_img_pt;
er_f11 = homog_img_pt;
homog_img_eid_diff = homog_img_pt;
homog_img_eid_abs = homog_img_pt;

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
extra={'ball1',[sb1]};

fmdl_t= ng_mk_cyl_models(0.,[16],[0.2,0,0.05],extra);
fmdl_t=fix_model(fmdl_t);
fmdl_t.stimulation = stim; %Add to model
fmdl_t.approx_type='tri3';
sim_img= mk_image(fmdl_t, cond_bkg );
pixel_group = fmdl_t.mat_idx{2};
%   figure; show_fem(homog_img);
%         pixel_group=[fmdl_t.mat_idx{2}];
%pixel_group = [102,123,103,83,66,82]; %b2C
%pixel_group = [327,364,328,292,259,291]; %c2C
%pixel_group = 1:length(fmdl.elems(:,1));

data = fwd_solve(sim_img);
data = add_noise(100, data, data_hom);



%% calc some derivatives
%obj_pixel_group = [100,120,144,138,115,95];
cond_obj = 2;


% Homog solution in recon domain
homog_img.fwd_solve.get_all_meas = 1;
vc = fwd_solve(homog_img); u0 = vc.volt;
DU0 = calc_grad_potential(homog_img, u0);

% Dataset
sim_img.elem_data(pixel_group) = cond_obj;
v0 = fwd_solve(sim_img);

% Where in image to vary conductivity
cond_try = 0.5:0.1:4;
%pixel_group = [100,120,144,138,115,95];

% Initialise
Ha = zeros(size(DU0,1),length(cond_try));
Hpf = Ha;
Hpd = Ha;

Ca = Ha;
Cpf = Ha;
Cpd = Ha;

Da = Ha;
Dpf = Ha;
Dpd = Ha;

% Clear cache so calc_hess_obj doesnt get confused
eidors_cache('clear_all')

% plot frequency
pl_frq = 5;
homog_img_out = homog_img;

%%

pixel_group = [655, 606, 607, 658, 677, 679, 678];

for ii=1:length(cond_try)

    % New delta_v
    homog_img.elem_data(pixel_group) = cond_try(ii);
    vi = fwd_solve(homog_img);
    delta_v = vi.meas - v0.meas;

    % Calculate with adjoint-state
    [H,D,C] = calc_hessian_obj(imdl.fwd_model, homog_img, 1:size(Ha,1), delta_v);
    Ha(:,ii)=diag(H);
    Da(:,ii)=diag(D);
    Ca(:,ii)=diag(C);
    
    % Calculate with PT and freespace Neumann func
    [Hpf(:,ii), Dpf(:,ii), Cpf(:,ii)] = calc_phessian_obj(imdl.fwd_model, homog_img, DU0, delta_v, 'freespace');
    

    % Calculate with PT and analytic disc Neumann func
    [Hpd(:,ii), Dpd(:,ii), Cpd(:,ii)] = calc_phessian_obj(imdl.fwd_model, homog_img, DU0, delta_v, 'disc');


end

%% Plots within the pixel groups for all cond vals
plot_pixel = [606, 481, 745]; 

for ii=1:length(plot_pixel)
    
    figure;
    plot(cond_try, Ha(plot_pixel(ii),:));
    hold all
    plot(cond_try, Hpf(plot_pixel(ii),:));
    plot(cond_try, Hpd(plot_pixel(ii),:));
    
%     figure;
%     plot(cond_try, Da(plot_pixel(ii),:));
%     hold all
%     plot(cond_try, Dpf(plot_pixel(ii),:));
%     plot(cond_try, Dpd(plot_pixel(ii),:));
%     
%     figure;
%     plot(cond_try, Ca(plot_pixel(ii),:));
%     hold all
%     plot(cond_try, Cpf(plot_pixel(ii),:));
%     plot(cond_try, Cpd(plot_pixel(ii),:));

end


%% Plots of whole domain
% find caxis range
cmin = min([Ca(:); Cpf(:); Cpd(:)]);
cmax = max([Ca(:); Cpf(:); Cpd(:)]);
homog_img.calc_colours.clim = [2.5e-9];

for ii=[4, 9, 14]
    figure;
    subplot(2,3,1)
    homog_img.elem_data = Ca(:,ii);
    show_fem(homog_img,1);
    
    
    subplot(2,3,2)
    homog_img.elem_data = Cpf(:,ii);
    show_fem(homog_img,1);
    
    subplot(2,3,3)
    homog_img.elem_data = Cpd(:,ii);
    show_fem(homog_img,1);

    subplot(2,3,4:6)
    plot(Ca(:,ii));
    hold all
    plot(Cpf(:,ii));
    plot(Cpd(:,ii));
    hold off
%     ylim([cmin, cmax])
    
end

