function inv_mdl= mk_common_model( str, varargin )
% MK_COMMON_MODEL: make common EIT models
%
% Utility function to create common EIT FEM models,
% so that users do not need to re-write common code
%
% Usage: 
%   mk_common_model('a2c',16)   - 2D circ model (64 elems) with 16 elecs
%   mk_common_model('b2c',16)   - 2D circ model (256 elems)
%   mk_common_model('c2c',16)   - 2D circ model (576 elems)
%   mk_common_model('d2c',16)   - 2D circ model (1024 elems)
%
%   mk_common_model('b3r1',16)  - circular ring with 16 electrodes
%   mk_common_model('b3r2',16)  - two circular rings with 16 electrodes
%
%   mk_common_model('b3z',16)   - zigzag pattern electrodes
%
%   mk_common_model('n3r2',16)  - NP's 3D model with 2 ring electrodes

options = {'no_meas_current','no_rotate_meas'};
n_elec= 16; % default

if     strcmp( str, 'a2c')
    inv_mdl = mk_2c_model( n_elec, 4, options );
elseif strcmp( str, 'b2c')
    inv_mdl = mk_2c_model( n_elec, 8, options );
elseif strcmp( str, 'c2c')
    inv_mdl = mk_2c_model( n_elec, 12, options );
elseif strcmp( str, 'd2c')
    inv_mdl = mk_2c_model( n_elec, 16, options );
elseif strcmp( str, 'b3z')
    inv_mdl = mk_dz_model( n_elec, options );
elseif strcmp( str, 'n3r2')
    inv_mdl = mk_n3r2_model( n_elec, options );
else
    error('don`t know what to do with option=',str);
end
    
function inv2d= mk_2c_model( n_elec, n_circles, options )

    n_rings= 1;
    params= mk_circ_tank(n_circles, [], n_elec); 

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

function inv3d= mk_dz_model( n_elec, options )

    n_rings= 1;

    levels= [-.4:.2:.4]; e_levels= [2,4]; nr= 8;
    levels= [-.5:.1:.5]; e_levels= [4,8]; nr= 4;

    params= mk_circ_tank( nr, levels, { 'zigzag', n_elec, e_levels } );

    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);

    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'np_fwd_solve';
    params.system_mat= 'np_calc_system_mat';
    params.jacobian=   'np_calc_jacobian';
    params.misc.sym= '{n}';
    fm3d = eidors_obj('fwd_model', params);

    inv3d.name=  'EIT inverse: 3D';
    %inv3d.solve= 'np_inv_solve';
     inv3d.solve= 'aa_inv_conj_grad'; % faster and feasible with less memory
    inv3d.hyperparameter.value = 1e-4;
    inv3d.image_prior.func= 'laplace_image_prior';
    inv3d.reconst_type= 'difference';
    inv3d.fwd_model= fm3d;
    inv3d= eidors_obj('inv_model', inv3d);

function inv_mdl = mk_n3r2_model( n_elec, options );
   load( 'datareal.mat' );
   bdy= dubs3( simp );
   fmdl.nodes= vtx;
   fmdl.elems= simp;
   fmdl.boundary= dubs3( simp );

   fmdl.solve=      'np_fwd_solve';
   fmdl.jacobian=   'np_calc_jacobian';
   fmdl.system_mat= 'np_calc_system_mat';

   for i=1:length(zc)
       electrodes(i).z_contact= zc(i);
       electrodes(i).nodes=     unique( elec(i,:) );
   end

   fmdl.gnd_node=           gnd_ind;
   fmdl.electrode =         electrodes;
   fmdl.misc.sym =          sym;

   [I,Ib] = set_3d_currents(protocol, elec, ...
               fmdl.nodes, fmdl.gnd_node, no_pl);

   % get the measurement patterns, only indH is used in this model
   %   here we only want to get the meas pattern from 'get_3d_meas',
   %   not the voltages, so we enter zeros
   [jnk,jnk,indH,indV,jnk] = get_3d_meas( elec, ...
            fmdl.nodes, zeros(size(I)), Ib, no_pl );
   n_elec= size(elec,1);
   n_meas= size(indH,1) / size(Ib,2);
   for i=1:size(Ib,2)
      fmdl.stimulation(i).stimulation= 'mA';
      fmdl.stimulation(i).stim_pattern= Ib(:,i);

      idx= ( 1+ (i-1)*n_meas ):( i*n_meas );
      fmdl.stimulation(i).meas_pattern= sparse ( ...
              (1:n_meas)'*[1,1], ...
              indH( idx, : ), ...
              ones(n_meas,2)*[1,0;0,-1], ...
              n_meas, n_elec );
   end
   fmdl= eidors_obj('fwd_model', fmdl);

   inv_mdl.name=         'Nick Polydorides EIT inverse';
   inv_mdl.solve=       'np_inv_solve';
   inv_mdl.hyperparameter.value = 1e-8;
   inv_mdl.image_prior.func= 'np_calc_image_prior';
   inv_mdl.image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
   inv_mdl.reconst_type= 'difference';
   inv_mdl.fwd_model= fmdl;
   inv_mdl= eidors_obj('inv_model', inv_mdl);
