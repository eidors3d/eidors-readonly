% How to make simulation data using EIDORS3D

% (C) 2005 Nick Polydorides + Andy Adler. Licenced under the GPL Version 2
% $Id: demo_3d_simdata.m,v 1.17 2005-10-28 15:10:55 aadler Exp $

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
params.solve=      'np_fwd_solve';
params.system_mat= 'np_calc_system_mat';
params.jacobian=   'np_calc_jacobian';
params.misc.perm_sym=   '{n}';
mdl_3d = eidors_obj('fwd_model', params);


disp('STEP 1A: simultion 3D - homogeneous');
% create homogeneous image + simulate data
cond= ones( size(mdl_3d.elems,1) ,1);
homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', cond, ...
                     'fwd_model', mdl_3d );
homg_data=fwd_solve( homg_img);

disp('STEP 1B: simultion 3D - inhomogeneous');

% create inhomogeneous image + simulate data
inhv= [38,50,51,66,67,83];
for inhlev= (e_levels(1)-1)*3 + [-3:2];
     cond(256*inhlev+inhv) =2;
%    cond(64*inhlev+[5,9,10,17,18,26]) =2;
end

inh_img= eidors_obj('image', 'inhomogeneous image', ...
                     'elem_data', cond, ...
                     'fwd_model', mdl_3d );
inh_data=fwd_solve( inh_img);
 show_fem( inh_img, 1); pause

% Add 10% noise
sig= std( inh_data.meas - homg_data.meas );
inh_data.meas= inh_data.meas + 0.10 * sig* randn( size(inh_data.meas) );

% 
% Step 2: Reconstruction in 2D
% 
params= mk_circ_tank(8, [], n_elec);

params.stimulation= stimulation;
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
 inv2d.image_prior.func= 'laplace_image_prior';
%inv2d.image_prior.func= 'aa_calc_image_prior';
inv2d.reconst_type= 'difference';
inv2d.fwd_model= mdl_2d_2;
inv2d= eidors_obj('inv_model', inv2d);

img2= inv_solve( inv2d, inh_data, homg_data);
img2.name= '2D inverse solution';
show_fem(img2); pause;

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
params.solve=      'np_fwd_solve';
params.system_mat= 'np_calc_system_mat';
params.jacobian=   'np_calc_jacobian';
params.misc.perm_sym= '{n}';
fm3d = eidors_obj('fwd_model', params);

inv3d.name=  'EIT inverse: 3D';
%inv3d.solve= 'np_inv_solve';
 inv3d.solve= 'aa_inv_conj_grad'; % faster and feasible with less memory
inv3d.hyperparameter.value = 1e-4;
inv3d.image_prior.func= 'laplace_image_prior';
inv3d.reconst_type= 'difference';
inv3d.fwd_model= fm3d;
inv3d= eidors_obj('inv_model', inv3d);

 img3= inv_solve( inv3d, inh_data, homg_data);
 img3.name= '3D inverse solution';
 show_fem(img3);
%figure;image_levels(img3, [-.4:.2:.4])
