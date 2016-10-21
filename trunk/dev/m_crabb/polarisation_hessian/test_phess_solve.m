% Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';



cond_bkg = 1;
cond_obj = 2;


%Make an image and show image
img = mk_image(fmdl,cond_bkg);

%figure; show_fem(img,[0,1,3])
img.fwd_solve.get_all_meas=1;

img.fwd_model.M_tensor.a = ones(size(img.elem_data,1),1);
img.fwd_model.M_tensor.b = ones(size(img.elem_data,1),1);
img.fwd_model.M_tensor.rot = zeros(size(img.elem_data,1),1);

img0=img;


pixel_group = [1,2,3,4];
img.elem_data(pixel_group) = cond_obj;

data = fwd_solve(img);



% Inverse
imdl.hyperparameter.value = 0.003;

% Default options
opts = [];
opts.ptensor = 1;
[x_P, r_P] = inv_solve_ptensor_lbfgs(imdl, img0, data, opts);

% Re-calc U0
opts.update_U0 = 1;
[x_Pup, r_Pup] = inv_solve_ptensor_lbfgs(imdl, img0, data, opts);

% Include reg
opts.update_U0 = 0;
opts.use_hyper = 1;
[x_Phy, r_Phy] = inv_solve_ptensor_lbfgs(imdl, img0, data, opts);

% Re-scaling
opts.use_hyper = 0;
opts.rescale = 1;
[x_Pre, r_Pre] = inv_solve_ptensor_lbfgs(imdl, img0, data, opts);

% no ptensor
opts.ptensor = 0;
[x_I, r_I] = inv_solve_ptensor_lbfgs(imdl, img0, data, opts);