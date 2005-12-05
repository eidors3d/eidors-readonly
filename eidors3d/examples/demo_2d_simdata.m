% Example of using EIDORS to simulate 2D data and to
% solve it using various 2D solvers

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: demo_2d_simdata.m,v 1.11 2005-12-05 22:12:11 aadler Exp $

% 
% Step 1: Create simple 16 electrode 2D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
n_elec= 16;
n_rings= 1;
%options = {'no_meas_current','rotate_meas'};
 options = {'no_meas_current','no_rotate_meas'};
params= mk_circ_tank(12, [], n_elec ); 

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
mdl_2d = eidors_obj('fwd_model', params);
show_fem( mdl_2d ); pause;

% create homogeneous image + simulate data
mat= ones( size(mdl_2d.elems,1) ,1);
homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_2d );
homg_data=fwd_solve( homg_img);

% create inhomogeneous image + simulate data
mat([65,81,82,101,102,122])=2;
inh_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_2d );
inh_data=fwd_solve( inh_img);

% 
% Step 2: Create different model for reconstruction
% 
params= mk_circ_tank(8, [], n_elec ); 

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
mdl_2d_2 = eidors_obj('fwd_model', params);
show_fem( mdl_2d_2 ); pause;

% 
% Step 3: Create inverse model
% 
clear inv2d;
inv2d.name= 'EIT inverse';
%inv2d.solve=       'aa_inv_solve';
 inv2d.solve=       'np_inv_solve';
%inv2d.solve=       'aa_inv_total_var';
 inv2d.hyperparameter.value = 1e-5;
%inv2d.hyperparameter.func = 'aa_calc_noise_figure';
%inv2d.hyperparameter.noise_figure= 2;
%inv2d.hyperparameter.tgt_elems= 1:4;
%inv2d.RtR_prior.func= 'laplace_image_prior';
 inv2d.R_prior.func= 'ab_calc_tv_prior';
%inv2d.RtR_prior.func= 'aa_calc_image_prior';
inv2d.reconst_type= 'difference';
inv2d.fwd_model= mdl_2d_2;
inv2d.fwd_model.misc.perm_sym= '{y}';
inv2d= eidors_obj('inv_model', inv2d);

% 
% Step 3: Reconst and show image
% 
img= inv_solve( inv2d, inh_data, homg_data);
show_slices(img);

