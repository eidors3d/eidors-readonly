% Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);




%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

imdl.fwd_model = fmdl;

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

data = fwd_solve(sim_img);

%% Eidors in-built inverse abs solver
imdl.hyperparameter.value = 0.00001;
imdl.fwd_model=fmdl;
imdl.solve =  @inv_solve_core;
imdl.reconst_type = 'absolute';
imdl.inve_solve_core = [];
img_eid_abs = inv_solve(imdl, data);

figure; 
subplot(121); show_fem(sim_img,1)
subplot(122); show_fem(img_eid_abs,1)

data_rec=fwd_solve(img_eid_abs)
figure; plot(data.meas,'r'); hold on; plot(data_hom.meas,'b');  hold on; plot(data_rec.meas,'b')

%% Hybrid update
u0 = data_hom.volt;
DU0 = calc_grad_potential(homog_img,u0);
imdl.inv_solve_core.update_func = @(J,W,dv,de,hp,opt)hybrid_GN_ptensor_update(sim_img.fwd_model,sim_img, DU0, J,W,dv,de,hp,opt);
img_hyb_abs = inv_solve(imdl, data);

figure; 
subplot(121); show_fem(sim_img,1)
subplot(122); show_fem(img_hyb_abs,1)

data_rec=fwd_solve(img_hyb_abs)
figure; plot(data.meas,'r'); hold on; plot(data_hom.meas,'b');  hold on; plot(data_rec.meas,'b')