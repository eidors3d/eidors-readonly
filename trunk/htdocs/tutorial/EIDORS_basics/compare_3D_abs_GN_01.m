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
v_i=fwd_solve(img_i); v_i_n = add_noise( 20, v_i ); v_h=fwd_solve(img_h);

%Plot actual and simulated voltage and show slice through 3D image
%figure; hold on; plot(v_i.meas); plot(v_h.meas,'r'); hold off;
clf; axes('position',[0.2,0.2,0.6,0.6]);

show_3d_slices(img_i,0.6,0.3,0.3); view(-24,12);
img_i.calc_colours.cb_shrink_move = [.3,.8,.00];
eidors_colourbar(img_i);
print_convert compare_3D_abs_GN_01a.png

show_fem(img_i);
print_convert compare_3D_abs_GN_01b.png
