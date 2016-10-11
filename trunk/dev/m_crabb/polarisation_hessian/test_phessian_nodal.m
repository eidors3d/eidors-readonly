%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

cond_bkg = 1;
cond_lower = 1;
cond_upper = 5;

cond_vals = linspace(cond_lower, cond_upper,5);

%Make an image and show image
img = mk_image(fmdl,cond_bkg);
img_backgrd=img; %Background to calculate Jacobians etc
%figure; show_fem(img,[0,1,3])

%Solve forward model on cond=1, save coefficients u0
%[d0, u0] = fwd_solve_higher_order_ptensor(img_backgrd);
img_backgrd.fwd_solve.get_all_meas=1;
[d0] = fwd_solve_higher_order(img_backgrd);
u0 = d0.volt;

% Gradient of solution on homogeneous domain
DU0 = calc_grad_potential(img_backgrd, u0);

%Pixk group of pixels to perturb
pixel_group = [327,364,328,292,259,291];

Hii = zeros(size(fmdl.elems,1), length(cond_vals));
PHii = Hii;
Gii = Hii;
D2ii = Hii;
PGii = Gii;
PD2ii = D2ii;

% M-tensor as circle
fmdl.M_tensor.a = ones(size(Hii,1),1);
fmdl.M_tensor.b = ones(size(Hii,1),1);
fmdl.M_tensor.rot = ones(size(Hii,1),1);
% img.elem_data(pixel_group) = 2;
% figure; show_fem(img,[1,0,0]);

for ii=1:length(cond_vals)
    img.elem_data(pixel_group)=cond_vals(ii);
    
        
    % Perturbed field
    %NB if want interior voltage here need to add
    %img.fwd_solve.get_all_meas=1 and use fwd_solve_higher_order.m
    d_perturb = fwd_solve(img);
    delta_d = d0.meas - d_perturb.meas;
    
    % phess
    
    [PHii(:,ii), PGii(:,ii), PD2ii(:,ii)] = calc_phessian_obj_OLD(fmdl, img_backgrd, DU0, [], delta_d,1);
    
    [H,G,D2] = calc_hessian_obj(fmdl, img_backgrd,1:size(Hii,1),delta_d);
    
    Hii(:,ii) = diag(H);
    Gii(:,ii) = diag(G);
    D2ii(:,ii) = diag(D2);
    
%     tic;
%     %Hii(:,ii) = 
%   foo =   calc_hessian_obj(fmdl, img, 1:size(Hii,1),delta_d);    
%     toc
    
end