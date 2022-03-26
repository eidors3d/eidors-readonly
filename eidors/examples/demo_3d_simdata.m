% How to make simulation data using EIDORS3D

% (C) 2005 Nick Polydorides + Andy Adler. License: GPL version 2 or version 3
% $Id$

% STIMULATION PATTERN
n_elec= 16;
n_rings= 1;
%options = {'no_meas_current','rotate_meas'};
 options = {'no_meas_current','no_rotate_meas'};
stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
% 
% Example 1: Create 16 electrode 3D model
% 
disp('STEP 1: Model simultion 3D');

% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
levels= [-.5:.1:.5];
e_levels= [4, 8];
% params= mk_circ_tank( 8, levels, n_elec );
  params= mk_circ_tank( 8, levels, { 'zigzag', n_elec, e_levels } );
% params= mk_circ_tank(12, levels, { 'zigzag', n_elec, [3,5,7] , ...
%                                    'planes', n_elec, 2} );

params.stimulation= stimulation;
params.solve=      @fwd_solve_1st_order;
%params.solve=      @eidors_default;
params.system_mat= @system_mat_1st_order;
params.jacobian=   @jacobian_adjoint;
mdl_3d = eidors_obj('fwd_model', params);


disp('STEP 1A: simultion 3D - homogeneous');
% create homogeneous image + simulate data
homg_img= mk_image(mdl_3d, 1);
homg_data=fwd_solve( homg_img);

disp('STEP 1B: simultion 3D - inhomogeneous');

% create inhomogeneous image + simulate data
inh_img = homg_img;
inhv= [38,50,51,66,67,83];
for inhlev= (e_levels(1)-1)*3 + [-3:2];
    inh_img.elem_data(256*inhlev+inhv) =2;
end
inh_data=fwd_solve( inh_img);
subplot(221);show_fem( inh_img);

% Add noise SNR=20
sig= std( inh_data.meas - homg_data.meas );
inh_data = add_noise(20, inh_data, homg_data);

% 
% Step 2: Reconstruction in 2D
% 
params= mk_circ_tank(8, [], n_elec);

params.stimulation= stimulation;
params.solve=      'fwd_solve_1st_order';
params.system_mat= 'system_mat_1st_order';
params.jacobian=   'jacobian_adjoint';
mdl_2d_2 = eidors_obj('fwd_model', params);

inv2d.name= 'EIT inverse';
inv2d.solve=       'inv_solve_diff_GN_one_step';
%inv2d.hyperparameter.func = 'select_noise_figure';
%inv2d.hyperparameter.noise_figure= 2;
%inv2d.hyperparameter.tgt_elems= 1:4;
 inv2d.hyperparameter.value = 1e-1;
 %inv2d.RtR_prior= 'prior_TV';
%inv2d.R_prior = 'prior_TV';
 inv2d.RtR_prior= 'prior_gaussian_HPF';
inv2d.jacobian_bkgnd.value= 1;
inv2d.reconst_type= 'difference';
inv2d.fwd_model= mdl_2d_2;
inv2d= eidors_obj('inv_model', inv2d);

img2= inv_solve( inv2d, homg_data, inh_data);
img2.name= '2D inverse solution';
subplot(223);   show_slices(img2);

% 
% Step 2: Reconstruction in 3D (using np_2003 code) and point
%          electrode models with 'zigzag' electrodes
% 
disp('STEP 2: Reconstruction 3D');
clear inv3d;

 levels= [-.4:.2:.4];
 params= mk_circ_tank( 8, levels, { 'zigzag', n_elec, [2,4] } );
%params= mk_circ_tank( 8, levels, { 'zigzag', n_elec, e_levels } );
%params= mk_circ_tank( 4, levels, { 'zigzag', n_elec, e_levels } );
%params= mk_circ_tank( 4, levels, n_elec );
params.stimulation= stimulation;
params.solve=      @fwd_solve_1st_order;
params.system_mat= @system_mat_1st_order;
params.jacobian=   @jacobian_adjoint;
params.misc.perm_sym= '{n}';
fm3d = eidors_obj('fwd_model', params);

inv3d.name=  'EIT inverse: 3D';
inv3d.solve= @inv_solve_diff_GN_one_step;
inv3d.hyperparameter.value = 1e-2;
inv3d.jacobian_bkgnd.value= 1;
%inv3d.RtR_prior= 'prior_TV';
inv3d.R_prior = 'prior_TV';
inv3d.reconst_type= 'difference';
inv3d.fwd_model= fm3d;
inv3d= eidors_obj('inv_model', inv3d);

 img3= inv_solve( inv3d, homg_data, inh_data);
 img3.name= '3D inverse solution';

subplot(122)
level(:,3) = [-.35:.2:.4]'; level(:,1:2) = Inf;
show_slices(img3, level);

