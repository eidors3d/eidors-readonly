n_elecs = 16;
if 0
   fmdl = ng_mk_cyl_models([2 2 0.2] ,[n_elecs,1],[0.1]);
else
   fmdl= mk_library_model('adult_male_16el');
   fmdl.electrode = fmdl.electrode([9:16,1:8]);
   fmdl.electrode = fmdl.electrode([1,16:-1:2]);
   fmdl.nodes(:,2) = -fmdl.nodes(:,2);
end
   fmdl.stimulation =  mk_stim_patterns(n_elecs,1,[0,1],[0,1],{'rotate_meas','no_meas_current'}, 1);

   fmdl.normalize_measurements = 1;
   img = mk_image(fmdl,1); % Homogeneous background

   opt.imgsz = [32 32];
   opt.distr = 3; 
   opt.noise_figure = 0.5; 
   imdl = mk_GREIT_model(img, 0.2, [], opt);

% MODEL
fmdl = imdl.rec_model;
fmdl.nodes = fmdl.nodes/max(fmdl.nodes(:))*1.2*[1,0;0,-1];
cmdl = mk_common_gridmdl('b2d','backproj');
cmdl = rmfield(cmdl.fwd_model,'coarse2fine');
f2c  = mk_coarse_fine_mapping(cmdl, fmdl);

RM = f2c*imdl.solve_use_matrix.RM(1:size(f2c,2),:);

ReconstrMatrix= - ( RM(1:2:end,:) + RM(2:2:end,:) )';
save ReconstrMatrixGREIT ReconstrMatrix -V6;
