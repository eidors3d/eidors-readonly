function inv_mdl= mk_common_model( str, n_elec, varargin )
% MK_COMMON_MODEL: make common EIT models
%
% Utility function to create common EIT FEM models,
% so that users do not need to re-write common code
%
% Usage: 
%      inv_mdl = mk_common_model( mdl_string, [n_elec/plane, n_planes])
%
% 2D Models circular models:
%   mk_common_model('a2c',16)   - 2D circ model (64 elems) with 16 elecs
%   mk_common_model('b2c',16)   - 2D circ model (256 elems)
%   mk_common_model('c2c',16)   - 2D circ model (576 elems)
%   mk_common_model('d2c',16)   - 2D circ model (1024 elems)
%   mk_common_model('e2c',16)   - 2D circ model (1600 elems)
%   mk_common_model('f2c',16)   - 2D circ model (2304 elems)
%
%   models ??c or ??c0 are rotated by zero.
%   models ??c1, ??c2, ??c3 are rotated by 22.5, 45, 67.5 degrees
%
% 2D Models using distmesh (D = show distmesh graphics, d= no graphics)
%   mk_common_model('a2d0c',16)  - 2D circ model using distmesh 
%   mk_common_model('b2d0c',16)  - 2D circ model using distmesh ~ 1300 elems
%   mk_common_model('d2d0c',16)  - 2D circ model using distmesh ~ 3200 elems
%
%      2d0 -> no 
%
% 2D Thorax models (levels 1 - 5 from shoulders to abdomen)
%   mk_common_model('b2t2',16)  - 2D Thorax#2 (chest) (256 elems)
%   mk_common_model('c2t4',16)  - 2D Thorax#3 (upper abdomen) (576 elems)
%   - all t1-t5 are available for each a-f models
%
% 2D square models:
%   mk_common_model('a2s',8)   - 2D square model (4x4x2 elems) (max 8 elecs)
%   mk_common_model('b2s',16)  - 2D square model (8x8x2 elems) (16 elecs)
%   mk_common_model('c2s',16)  - 2D square model (16x16x2 elems)
%   mk_common_model('d2s',16)  - 2D square model (24x24x2 elems)
%   mk_common_model('e2s',16)  - 2D square model (32x32x2 elems)
%   mk_common_model('f2s',16)  - 2D square model (40x40x2 elems)
%
%   models ??c or ??c0 are rotated by zero.
%   models ??c1, ??c2, ??c3 are rotated by 22.5, 45, 67.5 degrees
%
% 3D Models:
%   mk_common_model('n3r2',16)  - NP's 3D model with 2 ring electrodes
%   mk_common_model('n3z',16)   - NP's 3D model with zigzag electrodes
%
%   mk_common_model('b3cr',[16,3])  - cylinder with 3 rings of 16 elecs
%   mk_common_model('b3t2r',[16,1]) - t2 thorax shape with 1 ring of 16 elecs
%   mk_common_model('b3cz2',[16,1]) - cylinder with 2 rows of 8
%           zigzag pattern elecs. Stimulation treats this as 16x1 pattern
%   mk_common_model('b3cp2',16)      - cylinder with 2 rows of 8
%           elecs in 'planar' pattern. Stim treats this as 16x1 pattern
%
%   mk_common_model('a3cr',16)      - 64 elems * 4 planes
%   mk_common_model('b3cr',16)      - 256 elems * 10 planes 
%   mk_common_model('c3cr',16)      - 576 elems * 20 planes
%   mk_common_model('d3cr',16)      - 1024 elems * 40 planes
%

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_common_model.m,v 1.20 2008-03-29 00:42:11 aadler Exp $

options = {'no_meas_current','no_rotate_meas'};
% n_elec is number of [elec/ring n_rings]
if nargin<2
    n_elec= [16,1]; % default
end
if length(n_elec)==1
    n_elec= [n_elec,1]; % default
end
    
if length(str)<3
   error('format specified not recognized')
end

if str(2:3)=='2c'
% 2D circular models
   if     str(1)=='a'; layers=  4;
   elseif str(1)=='b'; layers=  8;
   elseif str(1)=='c'; layers= 12;
   elseif str(1)=='d'; layers= 16;
   elseif str(1)=='e'; layers= 20;
   elseif str(1)=='f'; layers= 24;
   else;  error(['don`t know what to do with option=%s',str]);
   end

   inv_mdl = mk_2c_model( n_elec, layers, options );   

   if length(str)==3; str= [str,'0'];end
      
   inv_mdl = rotate_model( inv_mdl, str2num(str(4)));

elseif str(2:3)=='2D'
   global distmesh_do_graphics; distmesh_do_graphics= 1;
   str= lower(str);
   inv_mdl= distmesh_2d_model(str, n_elec, options);
elseif str(2:3)=='2d'
   global distmesh_do_graphics; distmesh_do_graphics= 0;
   inv_mdl= distmesh_2d_model(str, n_elec, options);
elseif str(2:3)=='2s'
% 2D square models
   if     str(1)=='a'; layers=  4;
   elseif str(1)=='b'; layers=  8;
   elseif str(1)=='c'; layers= 16;
   elseif str(1)=='d'; layers= 24;
   elseif str(1)=='e'; layers= 32;
   elseif str(1)=='f'; layers= 40;
   else;  error('don`t know what to do with option=%s',str);
   end

   if rem( layers, n_elec(1)/2)~=0; 
      error('the %s model can`t support %d electrodes',str,n_elec(1));
   end
   inv_mdl = mk_2r_model( n_elec, layers, options);

elseif str(2:3)=='2t' & length(str)==4
   if     str(1)=='a'; layers=  4;
   elseif str(1)=='b'; layers=  8;
   elseif str(1)=='c'; layers= 12;
   elseif str(1)=='d'; layers= 16;
   elseif str(1)=='e'; layers= 20;
   elseif str(1)=='f'; layers= 24;
   else;  error(['don`t know what to do with option(1)=',str]);
   end

   inv_mdl = mk_2c_model( n_elec, layers, options );   
   inv_mdl = rotate_model( inv_mdl, 2); % 45 degrees

   if length(str)==0; str= [str,' '];end
      
   inv_mdl = deform_cylinder( inv_mdl, str2num(str(4)), 1 );

elseif str(2:3)=='3c' & length(str)==4
   if     str(1)=='a'; xy_layers=  4; z_layers= linspace(-.5,.5,5);
   elseif str(1)=='b'; xy_layers=  8; z_layers= linspace(-.7,.7,11);
   elseif str(1)=='c'; xy_layers= 12; z_layers= linspace(-.9,.9,21);
   elseif str(1)=='d'; xy_layers= 16; z_layers= linspace(-1,1,41);
   else;  error(['don`t know what to do with option(1)=',str]);
   end

   elec_per_plane = n_elec(1);
   spacing=.5;
   if     str(4)=='r';
      elec_conf= 'planes';
      ne_1  = (n_elec(2)-1)/2;
      elec_space= [ne_1:-1:-ne_1]*spacing;
   elseif str(4:5)=='z2';
      elec_conf= 'zigzag';
      elec_space= [1,-1]*spacing/2;
   elseif str(4:5)=='p2';
      elec_conf= 'planes';
      elec_space= [1,-1]*spacing/2;
      elec_per_plane = n_elec(1)/2;
   else;
      error(['don`t know what to do with option(4)=',str]);
   end

   inv_mdl = mk_3c_model( n_elec, xy_layers, z_layers, elec_space, ...
                              elec_per_plane, elec_conf, options );
%  inv_mdl = rotate_model( inv_mdl, n_elec(1)/8 );

elseif strcmp( str, 'n3r2')
    inv_mdl = mk_n3r2_model( n_elec, options );
elseif strcmp( str, 'n3z')
    inv_mdl = mk_n3z_model( n_elec, options );
else
    error(['Don`t know what to do with option=',str]);
end

inv_mdl.name= ['EIDORS common_model_',str]; 
inv_mdl= eidors_obj('inv_model', inv_mdl);
    
function inv2d = distmesh_2d_model(str, n_elec, options);
   if     str(1)=='a'; n_nodes=  50;
   elseif str(1)=='b'; n_nodes= 100;
   elseif str(1)=='c'; n_nodes= 200;
   elseif str(1)=='d'; n_nodes= 400;
   elseif str(1)=='e'; n_nodes= 800;
   elseif str(1)=='f'; n_nodes=1600;
   else;  error('don`t know what to do with option=%s',str);
   end

   elec_width= .1; refine_level=0;
   th=linspace(0,2*pi,n_elec(1)+1)';th(end)=[];
   elec_posn= [sin(th),cos(th)];
   [elec_nodes, refine_nodes] = dm_mk_elec_nodes( elec_posn, ...
          elec_width, refine_level);
   fd=inline('sqrt(sum(p.^2,2))-1','p');
   bbox = [-1,-1;1,1];
   z_contact = 0.01;
   fwd_mdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, ...
                             elec_nodes, refine_nodes, z_contact);

   inv2d= add_params_2d_mdl( fwd_mdl, n_elec(1), options);



function inv2d= mk_2c_model( n_elec, n_circles, options )

    n_elec= n_elec(1);
    params= mk_circ_tank(n_circles, [], n_elec); 
    inv2d= add_params_2d_mdl( params, n_elec, options);
   

function inv2d= mk_2r_model( n_elec, xy_size, options)
    if length(xy_size)==1; xy_size= xy_size*[1,1]; end
    xy_size= xy_size+1;

% TODO: To keep elements square, we should scale the 1's
    [x,y]= meshgrid( linspace(-1,1,xy_size(1)), ...
                     linspace(-1,1,xy_size(2)));
    x=x';y=y';
    fmdl.nodes= [x(:),y(:)];
    k= 1:xy_size(1)-1;
    elem_frac = [ k;k+1;k+xy_size(2); ...
                  k+1;k+xy_size(2);k+xy_size(2)+1];
    elem_frac= reshape(elem_frac, 3,[])';
    fmdl.elems=  [];
    for j=0:xy_size(2)-2
       fmdl.elems=  [fmdl.elems; elem_frac + xy_size(2)*j];
    end

    fmdl.boundary = find_boundary( fmdl.elems);

    % put 1/4 of elecs on each side 
    tb_elecs= linspace(1, xy_size(1), 1+2*n_elec(1)/4); 
    tb_elecs= tb_elecs(2:2:end);
    sd_elecs= linspace(1, xy_size(2), 1+2*n_elec(1)/4);
    sd_elecs= sd_elecs(2:2:end);
    
    el_nodes= [];
    % Top nodes -left to right
    bdy_nodes= (1:xy_size(1)) + xy_size(1)*(xy_size(2)-1); 
    el_nodes= [el_nodes, bdy_nodes(tb_elecs)];
    % Right nodes - top to bottom
    bdy_nodes= (1:xy_size(2))*xy_size(1); 
    el_nodes= [el_nodes, bdy_nodes(fliplr(sd_elecs))];
    % Bottom nodes - right to left
    bdy_nodes= 1:xy_size(1); 
    el_nodes= [el_nodes, bdy_nodes(fliplr(tb_elecs))];
    % Left nodes - bottom to top
    bdy_nodes= (0:xy_size(2)-1)*xy_size(1)+1; 
    el_nodes= [el_nodes, bdy_nodes(sd_elecs)];

%   trimesh(fmdl.elems,fmdl.nodes(:,1), fmdl.nodes(:,2));
    for i=1:n_elec(1)
       n= el_nodes(i);
       fmdl.electrode(i).nodes= n;
       fmdl.electrode(i).z_contact= .001; % choose a low value
%      plot(fmdl.nodes(n,1),fmdl.nodes(n,2),'*'); pause;
    end
    inv2d= add_params_2d_mdl( fmdl, n_elec(1), options);

% params is the part of the fwd_model
function inv2d= add_params_2d_mdl( params, n_elec, options);
    n_rings= 1;
    [st, els]= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', options, 10);
    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 0;
    params.np_fwd_solve.perm_sym= '{n}';
    mdl_2d   = eidors_obj('fwd_model', params);

    inv2d.solve=       'aa_inv_solve';
    inv2d.hyperparameter.value = 3e-2;
    %inv2d.hyperparameter.func = 'choose_noise_figure';
    %inv2d.hyperparameter.noise_figure= 1;
    %inv2d.hyperparameter.tgt_elems= 1:4;
     inv2d.RtR_prior= 'laplace_image_prior';
    %inv2d.RtR_prior= 'gaussian_HPF_prior';
    inv2d.jacobian_bkgnd.value= 1;
    inv2d.reconst_type= 'difference';
    inv2d.fwd_model= mdl_2d;

function inv3d = mk_3c_model( n_elec, xy_layers, z_layers, elec_space, ...
                              elec_per_plane, elec_conf, options );

    e_layers=[];
    for es= elec_space;
       ff= abs(z_layers  -es);
       ff= find(ff==min(ff));
       e_layers= [e_layers,ff(1)];
    end

    params= mk_circ_tank( xy_layers, z_layers, ...
           { elec_conf, elec_per_plane, e_layers} );

    [st, els]= mk_stim_patterns(n_elec(1), n_elec(2), '{ad}','{ad}', options, 10);

    params.stimulation= st;
    params.meas_select= els;
    params.solve=      'aa_fwd_solve';
    params.system_mat= 'aa_calc_system_mat';
    params.jacobian=   'aa_calc_jacobian';
    params.normalize_measurements= 0;
    params.np_fwd_solve.perm_sym= '{n}';
    fm3d = eidors_obj('fwd_model', params);

    inv3d.name=  'EIT inverse: 3D';
    inv3d.solve= @time_prior_solve;
    inv3d.time_prior_solve.time_steps=   0;
    inv3d.time_smooth_prior.space_prior = @noser_image_prior;
    inv3d.time_smooth_prior.time_weight = 0;
    inv3d.time_prior_solve.time_steps   = 0;

    inv3d.hyperparameter.value = 3e-2;
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
   fmdl.nodes= vtx;
   fmdl.elems= simp;
   fmdl.boundary= find_boundary( simp );

   fmdl.solve=      'np_fwd_solve';
   fmdl.jacobian=   'np_calc_jacobian';
   fmdl.system_mat= 'np_calc_system_mat';

   for i=1:length(zc)
       electrodes(i).z_contact= zc(i);
       electrodes(i).nodes=     unique( elec(i,:) );
   end

   fmdl.gnd_node=           gnd_ind;
   fmdl.electrode =         electrodes;
   fmdl.np_fwd_solve.perm_sym =     '{n}';

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
   inv_mdl.hyperparameter.value = 1e-2;
   inv_mdl.RtR_prior= 'np_calc_image_prior';
   inv_mdl.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
   inv_mdl.reconst_type= 'difference';
   inv_mdl.jacobian_bkgnd.value= 1;
   inv_mdl.fwd_model= fmdl;

% Deform the boundary of the cylinder to make it like a torso
% niv= 1.. 5 => Torso shape from T5 - T12
% xyz_expand - rescale xyz - default should be [1];
function inv_mdl = deform_cylinder( inv_mdl, niv, xyz_expand );
    NODE= inv_mdl.fwd_model.nodes';
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
      nn(1:2,i)= NODE(1:2,i)* fac;
    end  %for i=1:size

    inv_mdl.fwd_model.nodes = nn'*eye(xyz_expand);


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
    params.np_fwd_solve.perm_sym= '{n}';
    mdl_3d = eidors_obj('fwd_model', params);

    inv3d.name = 'EIT inverse: 3D';
    inv3d.solve=       'aa_inv_solve';
    %inv3d.solve=       'aa_inv_conj_grad';
    inv3d.hyperparameter.value = 1e-5;
    inv3d.RtR_prior= 'laplace_image_prior';
    %inv3d.RtR_prior= 'gaussian_HPF_prior';
    inv3d.jacobian_bkgnd.value= 1;
    inv3d.reconst_type= 'difference';
    inv3d.fwd_model= mdl_3d;

function inv3d= mk_b3r2_model( n_elec, nr, options )
    n_rings= 2;
    z_axis = [0:.1:1];
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
    params.np_fwd_solve.perm_sym= '{n}';
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
    %inv3d.RtR_prior= 'gaussian_HPF_prior';
    inv3d.jacobian_bkgnd.value= 1;
    inv3d.reconst_type= 'difference';
    inv3d.fwd_model= mdl_3d;

% rotate model rotate_model/16 times around
function inv_mdl = rotate_model( inv_mdl, rotate_mdl);
    nodes= inv_mdl.fwd_model.nodes;
    cos_rot = cos( pi/8*rotate_mdl);
    sin_rot = sin( pi/8*rotate_mdl);
    nodes(:,1)= inv_mdl.fwd_model.nodes(:,1:2)*[ cos_rot;-sin_rot];
    nodes(:,2)= inv_mdl.fwd_model.nodes(:,1:2)*[ sin_rot; cos_rot];
    inv_mdl.fwd_model.nodes= nodes;

    n_elec= length( inv_mdl.fwd_model.electrode );
    renum = rem(  rotate_mdl*n_elec/16 + (0:n_elec-1),n_elec)+1;
    inv_mdl.fwd_model.electrode = ...
       inv_mdl.fwd_model.electrode( renum);

