% Example of using EIDORS to simulate 2D data and to
% solve it using various 2D solvers

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

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
params.solve=      'fwd_solve_1st_order';
params.system_mat= 'system_mat_1st_order';
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
params.solve=      'fwd_solve_1st_order';
params.system_mat= 'system_mat_1st_order';
params.jacobian=   'jacobian_adjoint';
mdl_2d_2 = eidors_obj('fwd_model', params);
show_fem( mdl_2d_2 ); pause;

% 
% Step 3: Create inverse model
% 
clear inv2d;
inv2d.name= 'EIT inverse';
%inv2d.solve=       'inv_solve_diff_GN_one_step';
 inv2d.solve=       'np_inv_solve';
%inv2d.solve=       'aa_inv_total_var';
 inv2d.hyperparameter.value = 3e-3;
%inv2d.hyperparameter.func = 'select_noise_figure';
%inv2d.hyperparameter.noise_figure= 2;
%inv2d.hyperparameter.tgt_elems= 1:4;
%inv2d.RtR_prior= 'prior_laplace';
 inv2d.R_prior= 'prior_TV';
%inv2d.RtR_prior= 'prior_gaussian_HPF';
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;
inv2d.fwd_model= mdl_2d_2;
inv2d.fwd_model.misc.perm_sym= '{y}';
inv2d= eidors_obj('inv_model', inv2d);

% 
% Step 3: Reconst and show image
% 
img= inv_solve( inv2d, inh_data, homg_data);
show_slices(img);

