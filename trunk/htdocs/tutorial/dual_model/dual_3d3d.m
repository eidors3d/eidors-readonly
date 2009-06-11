% Dual model 3d-3d image reconstruction demo
%
% (C) 2009, Bartosz Sawicki
% Licenced under the GPLv2 or later
% $Id$

%% Create forward, fine tank model
electrodes_per_plane = 16;
number_of_planes = 2;
tank_radius = 0.2;
tank_height = 0.5;

[fine_mdl,centres] = create_tank_mesh_ng( ...
   tank_radius, ...   
   tank_height, ...  
   'C', ... % Rect_or_Circ_electrode
   log2(electrodes_per_plane), ...  
   number_of_planes, ...  
   0.15, ... % first_plane_starts
   0.2, ... % height_between_centres
   0.01, ... % electrode_width
   0.05, ... % electrode_height, ...
   'tank', ... % fname, 
   0 ... % elec_mesh_density  
   );

% Determine stimulation paterns
stim_pat = mk_stim_patterns(electrodes_per_plane, number_of_planes, ...
              '{ad}','{ad}',{'meas_current'});

% Parameters for forward model
fine_mdl.stimulation= stim_pat;
fine_mdl.solve=      'aa_fwd_solve';
fine_mdl.system_mat= 'aa_calc_system_mat';
fine_mdl.jacobian=   'aa_calc_jacobian';
fine_mdl.normalize_measurements= 0;
fine_mdl.np_fwd_solve.perm_sym= '{n}';

          
%% Solve homogeneous model 
disp('Homogeneous model');

% Every cells has the same material
mat= ones( size(fine_mdl.elems,1) ,1);

homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', fine_mdl );

homg_data=fwd_solve( fine_mdl, homg_img);


%% Create inclusion and solve inhomogenous model
disp('Inhomogeneous model');

% Parameters of spherical inclusion
center = [0.10, 0, 0.2];
radius = 0.05;
inclusion_material = 10;

inhomg_img = create_inclusion(homg_img, center, radius, inclusion_material);

figure; show_fem( inhomg_img );

inhomg_data=fwd_solve( fine_mdl, inhomg_img);


%% Create coarse model for inverse problem

coarse_mdl_maxh = 0.07; % maximum element size 
coarse_mdl = create_cylinder_mesh_ng('cylinder', tank_radius, ...
                                     tank_height, coarse_mdl_maxh) ;


%% Create inverse model
disp('Inverse model')

inv3d.name=  'Dual 3d-3d EIT inverse';
inv3d.fwd_model= fine_mdl;

disp('Calculating coarse2fine mapping ...');
inv3d.fwd_model.coarse2fine = ...
       mk_coarse_fine_mapping( fine_mdl, coarse_mdl);
disp('   ... done');

% Parameters for inverse model
inv3d.solve= 'np_inv_solve';
inv3d.hyperparameter.value = 1e-4;
inv3d.R_prior= 'np_calc_image_prior';
inv3d.np_calc_image_prior.parameters= [3 1]; 
inv3d.jacobian_bkgnd.value= 1;
inv3d.reconst_type= 'difference';

inv3d= eidors_obj('inv_model', inv3d);

recon_img= inv_solve( inv3d, homg_data, inhomg_data);

figure; show_fem( recon_img );


%% Show results on coarse mesh

mat= recon_img.elem_data;
coarse_recon_img= eidors_obj('image', 'reconstructed coarse image', ...
                     'elem_data', mat, ...
                     'fwd_model', coarse_mdl );

figure; show_fem( coarse_recon_img );              



