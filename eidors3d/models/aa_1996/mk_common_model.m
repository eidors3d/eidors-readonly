function inv_mdl= mk_common_model( str, varargin )
% MK_COMMON_MODEL: make common EIT models
%
% Utility function to create common EIT FEM models,
% so that users do not need to re-write common code
%
% Usage: 
%   mk_common_model('ac',16)   - 2D circ model with 16 electrodes
%
%   mk_common_model('dr',16)   - circular ring with 16 electrodes
%   mk_common_model('dr2',16)  - two circular rings with 16 electrodes
%
%   mk_common_model('dz',16)   - zigzag pattern electrodes

options = {'no_meas_current','no_rotate_meas'};
n_elec= 16; % default

if strcmp( str, 'ac')
    inv_mdl = mk_ac_model( n_elec, options );
else
    error('don`t know what to do with option=',str);
end
    
function inv2d= mk_ac_model( n_elec, options )

    n_rings= 1;
    params= mk_circ_tank(8, [], n_elec); 

    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 0;
    mdl_2d   = eidors_obj('fwd_model', params);

    inv2d.name= 'EIDORS model a0';
    inv2d.solve=       'aa_inv_solve';
    %inv2d.solve=       'aa_inv_conj_grad';
    inv2d.hyperparameter.value = 1e-5;
    %inv2d.hyperparameter.func = 'aa_calc_noise_figure';
    %inv2d.hyperparameter.noise_figure= 1;
    %inv2d.hyperparameter.tgt_elems= 1:4;
     inv2d.image_prior.func= 'laplace_image_prior';
    %inv2d.image_prior.func= 'aa_calc_image_prior';
    inv2d.reconst_type= 'difference';
    inv2d.fwd_model= mdl_2d;
    inv2d= eidors_obj('inv_model', inv2d);

