% 2D demo example for reconstruction of object floating inside tank with
% ambient fluid. Only internat structure of the object is reconstruted.
% External shape od object, and parameters of fluid is treated as a known.
% Dual meshing for inverse problem has been deployed.
% 
% (C) 2010 Bartosz Sawicki. License: GPL v.3

%% Create object model for reconstruction

object_center = [0.0, 0.1];
object_radius = 0.6;

mdl = create_circle_mesh_gmsh('circle', object_center, object_radius, 0.04 );
object_mdl = eidors_obj('fwd_model', mdl);

%% Create fine 'tank with object' model for forward problem

electrodes_per_plane = 16;
tank_radius = 1.0;

extra={'object','solid object = cylinder(0.0,0.,0;0,0,1;0.6) and plane(0,0,0;0,0,-1) and plane(0,0,0.2;0,0,1) -maxh=0.05;'};
[tank_mdl, mat_indices] = ng_mk_cyl_models([0,tank_radius,0.03],[electrodes_per_plane],[0.1,0,0.01], extra);

fluid_indices=mat_indices{1};
object_indices=mat_indices{2};

options = {'no_meas_current','no_rotate_meas'};
[st, els]= mk_stim_patterns(electrodes_per_plane, 1, '{ad}','{ad}', options, 10);

tank_mdl.stimulation= st;
tank_mdl.meas_select= els;
tank_mdl.solve=      'aa_fwd_solve';
tank_mdl.system_mat= 'bs_calc_system_mat';
tank_mdl.jacobian=   'aa_calc_jacobian';
tank_mdl.normalize_measurements= 0;
tank_mdl.np_fwd_solve.perm_sym= '{n}';

%% Solve model without inclusion
disp('Model without inclusion');

%Initialize
mat = ones(size(tank_mdl.elems,1),1);
%Fluid conductivity
mat(fluid_indices) = 0.75;
%Object avarage conductivity 
mat(object_indices) = 1.0;

homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', tank_mdl );

homg_data=fwd_solve( tank_mdl, homg_img);

%% Create inclusion and solve model
disp('Model with inclusion');

inclusion_center = [0.3, 0.3];
inclusion_radius = 0.05;
inclusion_conductivity = 1.1;

inhomg_img = create_inclusion(homg_img, inclusion_center, inclusion_radius, ...
    inclusion_conductivity);

inhomg_img.name = 'Inhomogeneous model';

figure; show_fem( inhomg_img );

inhomg_data = fwd_solve( tank_mdl, inhomg_img);

%% Add some noise
disp('Add noise');

SNR = -30; %dB

hdata = homg_data.meas;
idata = inhomg_data.meas;
noise = 10^(SNR/20)*std(idata-hdata)*randn(size(idata));
idata = idata + noise;
inhomg_data.meas = idata ;

%% Create inverse model
disp('Inverse model')

inv2d.name=  'Object in tank 2D EIT inverse';
inv2d.fwd_model= tank_mdl;
inv2d.rec_model= object_mdl;

disp('Calculating coarse2fine mapping ...');
inv2d.fwd_model.coarse2fine = ...
       mk_coarse_fine_mapping( tank_mdl, object_mdl);
disp('   ... done');

% Absolute reconstruction
inv2d.reconst_type= 'static';
inv2d.solve= @bs_nonlinearGN;
inv2d.parameters.term_tolerance= 1e-8;
inv2d.parameters.max_iterations= 5;
inv2d.nonlinearGN.init_backgnd= 0.00;
%
inv2d.RtR_prior= @tikhonov_image_prior;
inv2d.hyperparameter.value= 1e-3;
%
inv2d.fwd_model.background= homg_img.elem_data;

inv2d= eidors_obj('inv_model', inv2d);
recon_img= inv_solve(inv2d, inhomg_data);

figure; show_fem( recon_img );
