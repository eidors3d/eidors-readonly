function inv_mdl= mk_common_model( str, n_elec, varargin )
% MK_COMMON_MODEL: make common EIT models
%
% Utility function to create common EIT FEM models,
% so that users do not need to re-write common code
%
% Usage: 
% 3D Models:
%   mk_common_model('b3r1',16)  - circular ring with 16 electrodes
%   mk_common_model('b3r2',16)  - two circular rings with 16 electrodes
%
%   mk_common_model('b3z',16)   - zigzag pattern electrodes
%
%   mk_common_model('n3r2',16)  - NP's 3D model with 2 ring electrodes
%   mk_common_model('n3z',16)   - NP's 3D model with zigzag electrodes
%
% 2D Models:
%   mk_common_model('a2c',16)   - 2D circ model (64 elems) with 16 elecs
%   mk_common_model('b2c',16)   - 2D circ model (256 elems)
%   mk_common_model('c2c',16)   - 2D circ model (576 elems)
%   mk_common_model('d2c',16)   - 2D circ model (1024 elems)
%   mk_common_model('e2c',16)   - 2D circ model (1600 elems)
%   mk_common_model('f2c',16)   - 2D circ model (2304 elems)

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: mk_common_model.m,v 1.21 2006-08-28 12:29:32 aadler Exp $

options = {'no_meas_current','no_rotate_meas'};
% n_elec is number of [elec/ring n_rings]
if nargin<2
    n_elec= [16,1]; % default
end
    
if     strcmp( str, 'a2c')
    inv_mdl = mk_2c_model( n_elec, 4, options );
elseif strcmp( str, 'b2c')
    inv_mdl = mk_2c_model( n_elec, 8, options );
elseif strcmp( str, 'c2c')
    inv_mdl = mk_2c_model( n_elec, 12, options );
elseif strcmp( str, 'd2c')
    inv_mdl = mk_2c_model( n_elec, 16, options );
elseif strcmp( str, 'e2c')
    inv_mdl = mk_2c_model( n_elec, 20, options );
elseif strcmp( str, 'f2c')
    inv_mdl = mk_2c_model( n_elec, 24, options );   
elseif strcmp( str, 'b3z')
    inv_mdl = mk_dz_model( n_elec, options );
elseif strcmp( str, 'n3r2')
    inv_mdl = mk_n3r2_model( n_elec, options );
elseif strcmp( str, 'n3z')
    inv_mdl = mk_n3z_model( n_elec, options );
elseif strcmp( str, 'b3r1')
    inv_mdl = mk_b3r1_model( n_elec, options );
elseif strcmp( str, 'b3r2')
    inv_mdl = mk_b3r2_model( n_elec, options );    
else
    error('don`t know what to do with option=',str);
end

inv_mdl.name= ['EIDORS common_model_',str]; 
inv_mdl= eidors_obj('inv_model', inv_mdl);
    
function inv2d= mk_2c_model( n_elec, n_circles, options )

    n_elec= n_elec(1);
    n_rings= 1;
    params= mk_circ_tank(n_circles, [], n_elec); 

    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 0;
    params.misc.perm_sym= '{n}';
    mdl_2d   = eidors_obj('fwd_model', params);

    inv2d.solve=       'aa_inv_solve';
    %inv2d.solve=       'aa_inv_conj_grad';
    inv2d.hyperparameter.value = 3e-3;
    %inv2d.hyperparameter.func = 'aa_calc_noise_figure';
    %inv2d.hyperparameter.noise_figure= 1;
    %inv2d.hyperparameter.tgt_elems= 1:4;
     inv2d.RtR_prior= 'laplace_image_prior';
    %inv2d.RtR_prior= 'aa_calc_image_prior';
    inv2d.jacobian_bkgnd.value= 1;
    inv2d.reconst_type= 'difference';
    inv2d.fwd_model= mdl_2d;

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
    params.misc.perm_sym= '{n}';
    fm3d = eidors_obj('fwd_model', params);

    inv3d.name=  'EIT inverse: 3D';
    %inv3d.solve= 'np_inv_solve';
     inv3d.solve= 'aa_inv_conj_grad'; % faster and feasible with less memory
    inv3d.hyperparameter.value = 1e-4;
    inv3d.RtR_prior= 'laplace_image_prior';
    inv3d.reconst_type= 'difference';
    inv3d.jacobian_bkgnd.value= 1;
    inv3d.fwd_model= fm3d;

function inv_mdl = mk_n3z_model( n_elec, options );
   inv_mdl= mk_n3r2_model( n_elec, options);
   fwd_mdl= inv_mdl.fwd_model;
   renumber= [1:2:15; 18:2:32];
   fwd_mdl.electrode= fwd_mdl.electrode(renumber(:));
   n_rings= 1;
   [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
   fwd_mdl.stimulation= st;
   fwd_model.meas_select= els;
   inv_mdl.fwd_model= fwd_mdl;
   inv_mdl.name= 'NP 3D model with zigzag electrodes';

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
   fmdl.misc.perm_sym =     '{n}';

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
   inv_mdl.hyperparameter.value = 1e-4;
   inv_mdl.RtR_prior= 'np_calc_image_prior';
   inv_mdl.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
   inv_mdl.reconst_type= 'difference';
   inv_mdl.jacobian_bkgnd.value= 1;
   inv_mdl.fwd_model= fmdl;

% Deform the boundary of the cylinder to make it like a torso
% niv= 1.. 5 => Torso shape from T5 - T12
function inv_mdl = deform_cylinder( inv_mdl, niv, y_expand, x_expand)
    niv= reponse(2)-abs('0');
    x_coord= [ ...
      0,60,123,173,223,202,144,75,0,-75,-144,-202,-223,-173,-123,-60;
      0,50,105,138,144,144,109,50,0,-50,-109,-144,-144,-138,-105,-50;
      0,52, 99,133,148,141,110,61,0,-61,-110,-141,-148,-133,- 99,-52;
      0,51, 92,129,148,136, 96,47,0,-47,- 96,-136,-148,-129,- 92,-51;
      0,49, 92,128,148,141,111,64,0,-64,-111,-141,-148,-128,- 92,-49 ];
    y_coord= [ ...
      123,116, 91,42,-4,-67,-105,-119,-108,-119,-105,-67,-4,42, 91,116;
      129,132,112,62, 3,-57,-101,-110,-107,-110,-101,-57, 3,62,112,132;
      116,112, 92,53, 3,-48,- 88,-106,-105,-106,- 88,-48, 3,53, 92,112;
      143,130, 99,63,14,-35,- 68,- 82,- 82,- 82,- 68,-35,14,63, 99,130;
      136,128,103,68,23,-25,- 62,- 78,- 80,- 78,- 62,-25,23,68,103,128 ];

    reidx= [13:16, 1:12];
    geo= [x_coord(niv,reidx)',  ...
          y_coord(niv,reidx)'];
    a_max= size(geo,1);
    ab_geo=sqrt(sum(([ geo; geo(1,:) ]').^2)');
    nn= zeros(size(NODE));
    for i=1:size(NODE,2);
      angle = rem(a_max*atan2( NODE(2,i), ...
            NODE(1,i) )/2/pi+a_max,a_max)+1;
      fac=(  (floor(angle+1.001)- angle)* ...
              ab_geo(floor(angle+.001)) + ...
             (angle-floor(angle+.001))* ...
              ab_geo(floor(angle+1.001))  );
      nn(:,i)= NODE(:,i)* fac;
    end  %for i=1:size
    NODE=nn;

  NODE(1,:)= NODE(1,:)*x_expand;
  NODE(2,:)= NODE(2,:)*y_expand;
  NODE(3,:)= NODE(3,:)*z_expand;

function inv3d= mk_b3r1_model( n_elec, options )
    n_rings= 1;
    levels= [-.5:.1:.5]; 
    e_levels= 6; 
    nr= 8;

    params= mk_circ_tank( nr, levels, { 'planes', n_elec, e_levels } );
    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);

    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 0;
    params.misc.perm_sym= '{n}';
    mdl_3d = eidors_obj('fwd_model', params);

    inv3d.name = 'EIT inverse: 3D';
    inv3d.solve=       'aa_inv_solve';
    %inv3d.solve=       'aa_inv_conj_grad';
    inv3d.hyperparameter.value = 1e-5;
    inv3d.RtR_prior= 'laplace_image_prior';
    %inv3d.RtR_prior= 'aa_calc_image_prior';
    inv3d.jacobian_bkgnd.value= 1;
    inv3d.reconst_type= 'difference';
    inv3d.fwd_model= mdl_3d;

function inv3d= mk_b3r2_model( n_elec, options )
    n_rings= 2;
%    levels= [-.5:.1:.5]; 
    z_axis = [0:.1:1];  % show_slices() needs levels bw 0 and 1.
    e_levels= [4,8]; 
    nr= 4;
    n_elec = 8;
    
    params= mk_circ_tank( nr, z_axis, { 'planes', n_elec, e_levels } );
    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);

    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 0;
    params.misc.perm_sym= '{n}';
    mdl_3d = eidors_obj('fwd_model', params);
    
    % Specify number of levels in mesh for imaging slices
    num_levs = length(e_levels);
    levels = inf*ones(num_levs,3);
    levels(:,3) = e_levels / (length(z_axis)-1);
    levels(:,4) = ones(num_levs,1);
    levels(:,5) = (1:num_levs)';    
    mdl_3d.levels = levels;
    
    inv3d.name = 'EIT inverse: 3D';
    inv3d.solve=       'aa_inv_solve';
    %inv3d.solve=       'aa_inv_conj_grad';
    inv3d.hyperparameter.value = 1e-5;
    inv3d.RtR_prior= 'laplace_image_prior';
    %inv3d.RtR_prior= 'aa_calc_image_prior';
    inv3d.jacobian_bkgnd.value= 1;
    inv3d.reconst_type= 'difference';
    inv3d.fwd_model= mdl_3d;
