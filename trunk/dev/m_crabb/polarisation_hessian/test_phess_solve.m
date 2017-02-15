% Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);
imdl.fwd_model = fmdl;



%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';



cond_bkg = 1;
cond_obj = 2;


%Make an image and show image
sim_img = mk_image(fmdl,cond_bkg);

%figure; show_fem(img,[0,1,3])
sim_img.fwd_solve.get_all_meas=1;
sim_img.fwd_model.stimulation=stim;


sim_img.fwd_model.M_tensor.a = ones(size(sim_img.elem_data,1),1);
sim_img.fwd_model.M_tensor.b = ones(size(sim_img.elem_data,1),1);
sim_img.fwd_model.M_tensor.rot = zeros(size(sim_img.elem_data,1),1);

homog_img=sim_img;

data_hom = fwd_solve(homog_img);
%pixel_group = [1,2,3,4];
pixel_group = [115,138,95,137];
sim_img.elem_data(pixel_group) = cond_obj;

data = fwd_solve(sim_img)


%% Eidors in-built inverse diff solver
imdl.hyperparameter.value = 0.003;
img_eid_diff = inv_solve(imdl, data_hom,data);

figure; 
subplot(121); show_fem(sim_img,1)
subplot(122); show_fem(img_eid_diff,1)

%% Eidors in-built inverse abs solver
imdl.hyperparameter.value = 0.00001;
imdl.fwd_model=fmdl;
imdl.solve =  @inv_solve_core;
imdl.reconst_type = 'absolute';
img_eid_abs = inv_solve(imdl, data);

figure; 
subplot(121); show_fem(sim_img,1)
subplot(122); show_fem(img_eid_abs,1)

data_rec=fwd_solve(img_eid_abs)
figure; plot(data.meas,'r'); hold on; plot(data_hom.meas,'b');  hold on; plot(data_rec.meas,'b')

%%
% Default options
opts = [];
opts.ptensor = 1;
opts.use_hyper = 1;
opts.neumann = 'freespace';
[x_f, r_f] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);


%
opts.neumann = 'disc';
opts.use_hyper = 1;
[x_d, r_d] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);

% 
% % Re-calc U0
% opts.update_U0 = 1;
% [x_Pup, r_Pup] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% 
% % Include reg
% opts.update_U0 = 0;
% opts.use_hyper = 1;
% [x_Phy, r_Phy] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% 
% % Include reg and recalc U0
% opts.update_U0 = 1;
% opts.use_hyper = 1;
% [x_Phyup, r_Phyup] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% 
% % Re-scaling
% opts.use_hyper = 0;
% opts.rescale = 1;
% [x_Pre, r_Pre] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% 
% % no ptensor
% opts.ptensor = 0;
% [x_I, r_I] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);
% 



