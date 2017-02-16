% Calculate Hessian of objective, for a fixed dataset, as an element's
% contrast is varied.  See also comparess_obj_saturation, for varied
% dataset on fixed homogeneous domain.


% Data simulation model
smdl = mk_common_model('c2C',16);
fmdl = fix_model(smdl.fwd_model);
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stim = stim;
fmdl.approx_type = 'tri3';

smdl.fwd_model = fmdl;

% Inversion model
imdl = mk_common_model('b2C',16);
fmdl = fix_model(imdl.fwd_model);
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stim = stim;
fmdl.approx_type = 'tri3';


imdl.fwd_model = fmdl;


% imgs
cond_bkg = 1;

% Homog coarse img
img = mk_image(imdl, cond_bkg);
img.fwd_model.M_tensor.a = ones(size(img.elem_data,1),1);
img.fwd_model.M_tensor.b = ones(size(img.elem_data,1),1);
img.fwd_model.M_tensor.rot = zeros(size(img.elem_data,1),1);

homog_img = img;

% Fine img
sim_img = mk_image(smdl, cond_bkg);
sim_img.fwd_model.M_tensor.a = ones(size(sim_img.elem_data,1),1);
sim_img.fwd_model.M_tensor.b = ones(size(sim_img.elem_data,1),1);
sim_img.fwd_model.M_tensor.rot = zeros(size(sim_img.elem_data,1),1);

%% calc some derivatives
obj_pixel_group = [100,120,144,138,115,95];
cond_obj = 2;


% Homog solution in recon domain
img.fwd_solve.get_all_meas = 1;
vc = fwd_solve(img); u0 = vc.volt;
DU0 = calc_grad_potential(img, u0);

% Dataset
sim_img.elem_data(obj_pixel_group) = cond_obj;
v0 = fwd_solve(sim_img);

% Where in image to vary conductivity
cond_try = 0.5:0.1:4;
pixel_group = 45;

% Initialise
Ha = zeros(length(pixel_group),length(cond_try));
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

for ii=1:length(cond_try)

    % New delta_v
    img.elem_data(pixel_group) = cond_try(ii);
    vi = fwd_solve(img);
    delta_v = vi.meas - v0.meas;

    % Calculate with adjoint-state
    [HO, DO, CO] = calc_hessian_obj(imdl.fwd_model, img, pixel_group, delta_v);

    Ha(:,ii) = diag(HO);
    Ca(:,ii) = diag(CO);
    Da(:,ii) = diag(DO);
    
    % Calculate with PT and freespace Neumann func
    [HO, DO, CO] = calc_phessian_obj(imdl.fwd_model, img, DU0, delta_v, 'freespace');
    
    Hpf(:,ii) = HO(pixel_group);
    Dpf(:,ii) = DO(pixel_group);
    Cpf(:,ii) = CO(pixel_group);

    % Calculate with PT and analytic disc Neumann func
    [HO, DO, CO] = calc_phessian_obj(imdl.fwd_model, img, DU0, delta_v, 'disc');

    Hpd(:,ii) = HO(pixel_group);
    Dpd(:,ii) = DO(pixel_group);
    Cpd(:,ii) = CO(pixel_group);

end

figure;
plot(cond_try, Ca);
hold all
plot(cond_try, Cpf);
plot(cond_try, Cpd);

figure;
plot(cond_try, Da);
hold all
plot(cond_try, Dpf);
plot(cond_try, Dpd);

