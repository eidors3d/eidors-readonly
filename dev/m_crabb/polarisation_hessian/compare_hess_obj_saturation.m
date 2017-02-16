% Calculate Hessian of objective, at homogenous domain, with data from
% increasing contrast of object

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
sim_img.fwd_model.M_tensor.rot = zeros(sizeUntitled(sim_img.elem_data,1),1);


%% Loop through conductivities
cond_obj = 0.3:0.1:5;

pixel_group = [100,120,144,138,115,95];

% Initialise
Ha = zeros(size(img.elem_data,1),length(cond_obj));
Hpf = Ha;
Hpd = Ha;

Ca = Ha;
Cpf = Ha;
Cpd = Ha;

Da = Ha;
Dpf = Ha;
Dpd = Ha;


% 


% Homog data
v0 = fwd_solve(sim_img);

% Homog solution in reconstruction domain
img.fwd_solve.get_all_meas = 1;
vc = fwd_solve(img); u0 = vc.volt;
DU0 = calc_grad_potential(img, u0);

%%

for ii=1:length(cond_obj)
    
    sim_img.elem_data(pixel_group) = cond_obj(ii);
    vi = fwd_solve(sim_img);
    
    delta_v = vi.meas - v0.meas;
    
    % Calculate with adjoint-state
    [HO, DO, CO] = calc_hessian_obj(imdl.fwd_model, img, 1:size(img.elem_data,1), delta_v);
    
    Ha(:,ii) = diag(HO);
    Ca(:,ii) = diag(CO);
    Da(:,ii) = diag(DO);
    
    % Calculate with PT and freespace Neumann func
    [Hpf(:,ii), Dpf(:,ii), Cpf(:,ii)] = calc_phessian_obj(imdl.fwd_model, img, DU0, delta_v, 'freespace');
    
    % Calculate with PT and analytic disc Neumann func
    [Hpd(:,ii), Dpd(:,ii), Cpd(:,ii)] = calc_phessian_obj(imdl.fwd_model, img, DU0, delta_v, 'disc');
    
end
    
    
%% Now plot some results
% find caxis range
cmin = min([Ca(:); Cpf(:); Cpd(:)]);
cmax = max([Ca(:); Cpf(:); Cpd(:)]);
img.calc_colours.clim = [];

for ii=1:5:22
    figure;
    subplot(2,3,1)
    img.elem_data = Ca(:,ii);
    show_fem(img,1);
    
    
    subplot(2,3,2)
    img.elem_data = Cpd(:,ii);
    show_fem(img,1);
    
    subplot(2,3,3)
    img.elem_data = Cpd(:,ii);
    show_fem(img,1);

    subplot(2,3,4:6)
    plot(Ca(:,ii));
    hold all
    plot(Cpf(:,ii));
    plot(Cpd(:,ii));
    hold off
%     ylim([cmin, cmax])
    
end
    

%% Plot curve elem 95
elid = [161,114,75,44,21,11,42,110];%[144,100,64,36,16,1,27,85];

cmin = 1.1*min(min([Ca(elid,:); Cpf(elid,:); Cpd(elid,:)]));
cmax = 1.5*max(max([Ca(elid,:); Cpf(elid,:); Cpd(elid,:)]));

figure
for ii=1:8
    subplot(2,4,ii)
    plot(cond_obj,Ca(elid(ii), :));
    grid on
    hold all
    plot(cond_obj,Cpf(elid(ii),:));
    plot(cond_obj,Cpd(elid(ii),:));
    hold off
    
    ylim([cmin, cmax])

end


