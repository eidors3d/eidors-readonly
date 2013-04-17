%Homogeneous and inclusion conductivity
cond_h=1.0; cond_inc=2.0;

%3D forward model without inclusion 
mdl_h= ng_mk_cyl_models([1 .7],[16,.25,.75],[0.075,0.3]);
mdl_h.normalize_measurements=0;

%3D forward model with inclusion 
extra={'ball','solid ball = sphere(0.2,0.2,0.5;0.2);'};
[mdl_i,mat_idx_i]= ng_mk_cyl_models([1 .7],[16,.25,.75],[0.075,0.3],extra);
mdl_i.normalize_measurements=0;

%Stimulation patterns and add to models
stim=mk_stim_patterns(16,2,'{ad}','{ad}');
mdl_h.stimulation=stim; mdl_i.stimulation=stim;

%Create two images
img_h= mk_image(mdl_h,cond_h); 
img_i= mk_image(mdl_i,cond_h); img_i.elem_data(mat_idx_i{2}) = cond_inc;

%Now get "real" voltages and add noise
v_i=fwd_solve(img_i); v_i_n = add_noise( 10, v_i ); v_h=fwd_solve(img_h);

%Plot actual and simulated voltage and show slice through 3D image
%figure; hold on; plot(v_i.meas); plot(v_h.meas,'r'); hold off;
figure; show_3d_slices(img_i,0.6,0.3,0.3); eidors_colourbar(img_i);


%Inverse solution
%Create generic absolute reconstruction model
imdl = mk_common_model('b2c2',32);
imdl.fwd_model = mdl_h;
imdl.reconst_type = 'absolute';
imdl.jacobian_bkgnd.value=img_h.elem_data; %Background conductivity
imdl.parameters.show_iterations=2; %Show iteration progress
imdl.parameters.best_homog=2; %Best fitting homogeneous
imdl.parameters.max_iterations = 10; %Number of iterations



%(i) Tikhonov prior
%Loop over hyperparameter, prior types and iterations
imdl.RtR_prior=@prior_tikhonov;


%Default Gauss Newton solvers
imdl.hyperparameter.value = 0.003^2;
imdl.solve = @inv_solve_abs_GN;
img = inv_solve(imdl , v_i); img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img_n = inv_solve(imdl  , v_i_n); img_n.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
vr_gn=fwd_solve(img); vr_gn_n=fwd_solve(img_n);
figure; hold on; plot(v_i.meas-vr_gn.meas);  plot(v_i_n.meas-vr_gn_n.meas,'r');
figure; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
figure; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);

%Alternative Gauss Newton solver, changing prior at each iteration
imdl.hyperparameter.value = 0.003;
imdl.solve = @inv_solve_abs_GN_prior;
img = inv_solve(imdl , v_i); img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img_n = inv_solve(imdl  , v_i_n); img_n.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
vr_agn=fwd_solve(img); vr_agn_n=fwd_solve(img_n);
figure; hold on; plot(v_i.meas-vr_agn.meas);  plot(v_i.meas-vr_agn_n.meas,'g'); plot(v_i.meas,'r'); hold off;
figure; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
figure; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);

%Constrained Gauss Newton solver
imdl.parameters.min_s=0.8; %Minimum conductivity
imdl.parameters.max_s=2.2; %Maximum conductivity
imdl.parameters.rel_par=2.0; %Surrogate parameter
imdl.hyperparameter.value = 0.003;
imdl.solve = @inv_solve_abs_GN_constrain;
img = inv_solve(imdl , v_i); img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img_n = inv_solve(imdl  , v_i_n); img_n.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
vr_cgn=fwd_solve(img); vr_cgn_n=fwd_solve(img_n);
figure; hold on; plot(v_i.meas-vr_cgn.meas); plot(v_i.meas-vr_cgn_n.meas,'g'); plot(v_i.meas,'r'); hold off;
figure; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
figure; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);



%(ii) NOSER prior 
%Loop over hyperparameter, prior types and iterations
imdl.RtR_prior=@prior_noser;

%Default Gauss Newton solvers
imdl.hyperparameter.value = 0.03^2;
imdl.solve = @inv_solve_abs_GN;
img = inv_solve(imdl , v_i); img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img_n = inv_solve(imdl  , v_i_n); img_n.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
vr_gn=fwd_solve(img); vr_gn_n=fwd_solve(img_n);
figure; hold on; plot(v_i.meas-vr_gn.meas);  plot(v_i.meas-vr_gn_n.meas,'g'); plot(v_i.meas,'r'); hold off;
figure; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
figure; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);

%Alternative Gauss Newton solver, changing prior at each iteration
imdl.hyperparameter.value = 0.03;
imdl.solve = @inv_solve_abs_GN_prior;
img = inv_solve(imdl , v_i); img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img_n = inv_solve(imdl  , v_i_n); img_n.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
vr_agn=fwd_solve(img); vr_agn_n=fwd_solve(img_n);
figure; hold on; plot(v_i.meas-vr_agn.meas);  plot(v_i.meas-vr_agn_n.meas,'g'); plot(v_i.meas,'r'); hold off;
figure; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
figure; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);

%Constrained Gauss Newton solver
imdl.parameters.min_s=0.8; %Minimum conductivity
imdl.parameters.max_s=2.2; %Maximum conductivity
imdl.hyperparameter.value = 0.03;
imdl.solve = @inv_solve_abs_GN_constrain;
img = inv_solve(imdl , v_i); img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
img_n = inv_solve(imdl  , v_i_n); img_n.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
vr_cgn=fwd_solve(img); vr_cgn_n=fwd_solve(img_n);
figure; hold on; plot(v_i.meas-vr_cgn.meas); plot(v_i.meas-vr_cgn_n.meas,'g'); plot(v_i.meas,'r'); hold off;
figure; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
figure; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);
