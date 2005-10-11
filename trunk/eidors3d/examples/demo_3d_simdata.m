% How to make simulation data using EIDORS3D
% $Id: demo_3d_simdata.m,v 1.6 2005-10-11 20:27:04 aadler Exp $

% 
% Example 1: Create simple 16 electrode 2D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
n_elec= 16;
n_rings= 1;
levels= [-.5:.1:.5];
% params= mk_circ_tank(12, levels, n_elec );
  params= mk_circ_tank(12, levels, { 'zigzag', n_elec, [4,8] } );
% params= mk_circ_tank(12, levels, { 'zigzag', n_elec, [3,5,7] , ...
%                                    'planes', n_elec, 2} );

%options = {'no_meas_current','rotate_meas'};
 options = {'no_meas_current','no_rotate_meas'};
params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
mdl_3d = eidors_obj('fwd_model', params);
 show_fem( mdl_3d ); pause


% create homogeneous image + simulate data
mat= ones( size(mdl_3d.elems,1) ,1);
homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_3d );
homg_data=fwd_solve( homg_img);

% create inhomogeneous image + simulate data
mat(50)= 2;
inh_img= eidors_obj('image', 'inhomogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_3d );
inh_data=fwd_solve( inh_img);
figure;image_levels(inh_img, [-.4:.2:.4])

% 
% Step 2: Reconstruction in 3D (using np_2003 code)
% 
clear inv3d;
params= mk_circ_tank(4 , levels, { 'zigzag', n_elec, [4,8] } );
params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
params.misc.sym= '{n}';
fm3d = eidors_obj('fwd_model', params);

inv3d.name=  'EIT inverse: 3D';
inv3d.solve= 'np_inv_solve';
inv3d.hyperparameter.value = 1e-2;
inv3d.image_prior.func= 'tikhonov_image_prior';
inv3d.reconst_type= 'difference';
inv3d.fwd_model= fm3d;
inv3d= eidors_obj('inv_model', inv3d);

img3= inv_solve( inv3d, inh_data, homg_data);
%show_fem(img3);
figure;image_levels(img3, [-.4:.2:.4])

% 
% Step 3: Reconstruction in 2D
% 
params= mk_circ_tank(8, [], n_elec);

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
mdl_2d_2 = eidors_obj('fwd_model', params);

inv2d.name= 'EIT inverse';
inv2d.solve=       'aa_inv_solve';
%inv2d.hyperparameter.func = 'aa_calc_noise_figure';
%inv2d.hyperparameter.noise_figure= 2;
%inv2d.hyperparameter.tgt_elems= 1:4;
 inv2d.hyperparameter.value = 1e-2;
 inv2d.image_prior.func= 'tikhonov_image_prior';
%inv2d.image_prior.func= 'aa_calc_image_prior';
inv2d.reconst_type= 'difference';
inv2d.fwd_model= mdl_2d_2;
inv2d= eidors_obj('inv_model', inv2d);

img2= inv_solve( inv2d, inh_data, homg_data);
show_fem(img);

