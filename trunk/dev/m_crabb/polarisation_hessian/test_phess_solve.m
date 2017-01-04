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

sim_img.fwd_model.M_tensor.a = ones(size(sim_img.elem_data,1),1);
sim_img.fwd_model.M_tensor.b = ones(size(sim_img.elem_data,1),1);
sim_img.fwd_model.M_tensor.rot = zeros(size(sim_img.elem_data,1),1);

homog_img=sim_img;


pixel_group = [1,2,3,4];
sim_img.elem_data(pixel_group) = cond_obj;

data = fwd_solve(sim_img);



% Inverse
imdl.hyperparameter.value = 0.003;
imdl.reconst_type = 'static';

%% Eidors in-built
img_eid = inv_solve(imdl, data);


%% 

% Default options
opts = [];
opts.ptensor = 1;
[x_P, r_P] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);

% Re-calc U0
opts.update_U0 = 1;
[x_Pup, r_Pup] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);

% Include reg
opts.update_U0 = 0;
opts.use_hyper = 1;
[x_Phy, r_Phy] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);

% Re-scaling
opts.use_hyper = 0;
opts.rescale = 1;
[x_Pre, r_Pre] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);

% no ptensor
opts.ptensor = 0;
[x_I, r_I] = inv_solve_ptensor_lbfgs(imdl, homog_img, data, opts);




