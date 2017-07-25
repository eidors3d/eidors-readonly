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
pixel_group = [100,120,144,138,115,95];

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
img_out = img;

%%

for ii=1:length(cond_try)

    % New delta_v
    img.elem_data(pixel_group) = cond_try(ii);
    vi = fwd_solve(img);
    delta_v = vi.meas - v0.meas;

    % Calculate with adjoint-state
    [H,D,C] = calc_hessian_obj(imdl.fwd_model, img, 1:size(Ha,1), delta_v);
    Ha(:,ii)=diag(H);
    Da(:,ii)=diag(D);
    Ca(:,ii)=diag(C);
    
    % Calculate with PT and freespace Neumann func
    [Hpf(:,ii), Dpf(:,ii), Cpf(:,ii)] = calc_phessian_obj(imdl.fwd_model, img, DU0, delta_v, 'freespace');
    

    % Calculate with PT and analytic disc Neumann func
    [Hpd(:,ii), Dpd(:,ii), Cpd(:,ii)] = calc_phessian_obj(imdl.fwd_model, img, DU0, delta_v, 'disc');


end

%% Plots within the pixel groups for all cond vals
for ii=1:length(pixel_group)
    
    figure;
    plot(cond_try, Ca(pixel_group(ii),:));
    hold all
    plot(cond_try, Cpf(pixel_group(ii),:));
    plot(cond_try, Cpd(pixel_group(ii),:));
    
    figure;
    plot(cond_try, Da(pixel_group(ii),:));
    hold all
    plot(cond_try, Dpf(pixel_group(ii),:));
    plot(cond_try, Dpd(pixel_group(ii),:));

end


%% Plots of whole domain
% find caxis range
cmin = min([Ca(:); Cpf(:); Cpd(:)]);
cmax = max([Ca(:); Cpf(:); Cpd(:)]);
img.calc_colours.clim = [];

for ii=[4, 9, 19]
    figure;
    subplot(2,3,1)
    img.elem_data = Ca(:,ii);
    show_fem(img,1);
    
    
    subplot(2,3,2)
    img.elem_data = Cpf(:,ii);
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

