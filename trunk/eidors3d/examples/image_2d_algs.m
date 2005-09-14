% Based on the 'bubble' data from Eidors2D, use several 
% different algorithms to image it
%
% $Id: image_2d_algs.m,v 1.2 2005-09-14 22:15:12 aadler Exp $

eidors_msg('log_level',1); % 2 for most messages

% 
% Step 1: Create simple 16 electrode 2D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
n_elec= 16;
n_rings= 1;
 options = {'no_meas_current','no_rotate_meas'};

params= mk_circ_tank(8, [], n_elec, n_rings); 

[st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
params.stimulation= st;
params.meas_select= els;
params.solve=      'aa_fwd_solve';
params.system_mat= 'aa_calc_system_mat';
params.jacobian=   'aa_calc_jacobian';
params.normalize_measurements= 1;
mdl_2d   = eidors_obj('fwd_model', params);

% 
% Step 3: Create inverse model
% 
inv2d.name= 'AA mdl with excluded measurements';
inv2d.solve=       'aa_inv_solve';
 inv2d.hyperparameter.value = 1e+2;
%   inv2d.hyperparameter.func = 'aa_calc_noise_figure';
%   inv2d.hyperparameter.noise_figure= 1;
%   inv2d.hyperparameter.tgt_elems= 1:4;
 inv2d.image_prior.func= 'tikhonov_image_prior';
%inv2d.image_prior.func= 'aa_calc_image_prior';
inv2d.reconst_type= 'difference';
inv2d.fwd_model= mdl_2d;
inv2d= eidors_obj('inv_model', inv2d);

% 
% Step 3: Reconst and show image
% 
load eidors2d_bubble.mat
d1= bubble2(1280+(-255:0));
d2= bubble1(1280+(-255:0));

img= inv_solve( inv2d, d1, d2);
show_fem(img);

