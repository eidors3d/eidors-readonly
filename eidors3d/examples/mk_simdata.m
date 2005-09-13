% How to make simulation data using EIDORS3D
% $Id: mk_simdata.m,v 1.12 2005-09-13 00:56:58 aadler Exp $

eidors_msg('log_level',1); % 2 for most messages

% 
% Example 1: Create simple 16 electrode 2D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
n_elec= 16;
n_rings= 1;
%options = {'no_meas_current','rotate_meas'};
 options = {'no_meas_current','no_rotate_meas'};
params= mk_circ_tank(8, [], n_elec, n_rings); 

params.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            options, 10);
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
mdl_2d = eidors_obj('fwd_model', params);
show_fem( mdl_2d );

% create homogeneous image + simulate data
mat= ones( size(mdl_2d.elems,1) ,1);
homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_2d );
homg_data=fwd_solve( homg_img);

% create inhomogeneous image + simulate data
mat(5)= 2;
inh_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_2d );
inh_data=fwd_solve( inh_img);

inv2d.name= 'EIT inverse';
inv2d.solve=       'aa_inv_solve';
%inv2d.hyperparameter.value = 1e-8;
inv2d.hyperparameter.func = 'aa_calc_noise_figure';
inv2d.hyperparameter.noise_figure= 2;
inv2d.hyperparameter.tgt_elems= 1:4;
inv2d.image_prior.func= 'tikhonov_image_prior';
inv2d.reconst_type= 'difference';
inv2d.fwd_model= mdl_2d;
inv2d= eidors_obj('inv_model', inv2d);

img= inv_solve( inv2d, inh_data, homg_data);

% 
% Example 2: Create simple 16 electrode 3D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
params= mk_circ_tank(8, [-.5,0,.5], n_elec, n_rings); 
mdl_3d = eidors_obj('fwd_model', params);
