% STEP 4: Reconstruction using GREIT 3D
   vopt.imgsz = [32 32];
   vopt.zvec = linspace( 0.75,3.25,6);
   vopt.save_memory = 1;
   opt.noise_figure = 2.0;
   [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
   imdl3= mk_GREIT_model(imdl, 0.20, [], opt);
