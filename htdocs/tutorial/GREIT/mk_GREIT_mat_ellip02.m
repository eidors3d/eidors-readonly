% $Id$
opt.imgsz = [32 32];
opt.distr = 3; % non-random, uniform
opt.Nsim = 1000;
opt.target_size = 0.05; % Target size (frac of medium)
opt.noise_figure = 0.5; % Recommended NF=0.5;

k=1; for el = linspace(1,2.4,4);
   fmdl = ng_mk_ellip_models([2, 2,el,0.2] ,[n_elecs,1],[0.1]);
   fmdl.stimulation =  stim;
   fmdl.normalize_measurements = 0;
   img = mk_image(fmdl,1);

   imdl(k) = mk_GREIT_model(img, 0.25, [], opt);

   subplot(1,4,k); show_fem(fmdl); view(0,70);
   k=k+1;
end

print_convert mk_GREIT_mat_ellip02.png '-density 180'
