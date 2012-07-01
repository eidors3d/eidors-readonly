n_elecs = 16;
   fmdl= mk_library_model('adult_male_16el');
   fmdl.electrode = fmdl.electrode([9:16,1:8]);
   fmdl.electrode = fmdl.electrode([1,16:-1:2]);
   fmdl.stimulation =  mk_stim_patterns(n_elecs,1,[0,1],[0,1],{'rotate_meas','no_meas_current'}, 1);

   fmdl = mdl_normalize(fmdl, 1);
   img = mk_image(fmdl,1); % Homogeneous background

   opt.imgsz = [32 32];
   opt.distr = 3; 
   opt.noise_figure = 0.5; 
   imdl = mk_GREIT_model(fmdl, 0.25, [], opt);

% MODEL
fmdl = imdl.rec_model;
fmdl.nodes = fmdl.nodes/max(fmdl.nodes(:))*1.1*[1,0;0,-1];
cmdl = mk_common_gridmdl('b2d','backproj');
cmdl = rmfield(cmdl.fwd_model,'coarse2fine');
f2c  = mk_coarse_fine_mapping(cmdl, fmdl);

RM = f2c*imdl.solve_use_matrix.RM(1:size(f2c,2),:);

ReconstrMatrix= - ( RM(1:2:end,:) + RM(2:2:end,:) )';
save ReconstrMatrixGREITt ReconstrMatrix -V6;
