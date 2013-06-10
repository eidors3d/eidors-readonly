n_elecs = 16;
fmdl = ng_mk_cyl_models([1 1 0.1] ,[n_elecs,.5],[0.05]);
% NOTE: The 'rotate_meas' is essential here.
fmdl.stimulation =  mk_stim_patterns(n_elecs,1,[0,1],[0,1],{'rotate_meas','no_meas_current'}, 1);
fmdl = mdl_normalize(fmdl, 1); %MUST FOR AUSPEX

opt.imgsz = [32 32];
opt.distr = 3; 
opt.noise_figure = 0.5; 
imdl = mk_GREIT_model(fmdl, 0.20, [], opt);

% DESTINATION AUSPEX MODEL
cmdl = mk_common_gridmdl('b2d','backproj');
cmdl = rmfield(cmdl.fwd_model,'coarse2fine');
% MODEL GEOMETRIES MUST MATCH
f2c  = mk_coarse_fine_mapping(cmdl, imdl.rec_model);

RM = f2c*imdl.solve_use_matrix.RM(1:size(f2c,2),:);

ReconstrMatrix= - ( RM(1:2:end,:) + RM(2:2:end,:) )';
save ReconstrMatrixGREITc ReconstrMatrix -V6;
