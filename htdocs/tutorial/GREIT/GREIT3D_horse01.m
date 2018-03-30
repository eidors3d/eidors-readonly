skip5 = {32,1,[0,5],[0,5],{'no_meas_current_next1'},1};

fmdl= ng_mk_ellip_models([4,0.8,1.1,.5],[16,1.7,2.3],[0.05]);
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip5{:});

% From carleton/2017/conference/EIT2017/3D-EIT-Horses/code/data_reconst.m
  extraflip= [4:12];
  idx = reshape(1:32,2,[])';
  idx(2:2:end,:) = fliplr(idx(2:2:end,:));
  idx(extraflip,:) = fliplr(idx(extraflip,:));
  fmdl.electrode(idx) = fmdl.electrode(:);

   vopt.imgsz = [32 32];
   vopt.square_pixels = true;
   vopt.zvec = linspace(-1,1,10)*1.2+2;
   vopt.save_memory = 1;
   opt.noise_figure = 1.0;
   [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
   imdl3= mk_GREIT_model(imdl, 0.20, [], opt);


fmdl= ng_mk_ellip_models([4,0.8,1.1,.5],[32,2.0],[0.05]);
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip5{:});

   [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
   imdl2a= mk_GREIT_model(imdl, 0.20, [], opt);

   clear opt;
   opt.imgsz = [32 32];
   opt.square_pixels = true;
   opt.noise_figure = 0.5;
   img = mk_image(fmdl,1);
   imdl2b= mk_GREIT_model(img, 0.25, [], opt);

load mar28-horse_breathing.mat

for i=1:3; switch i;
   case 1; imdl = imdl3;  vv= horse3d; % 2 planes, 3D GREIT
   case 2; imdl = imdl2a; vv= horse2d; % 1 plane, 3D GREIT
   case 3; imdl = imdl2b; vv= horse2d; % 1 plane, 2D GREIT
   end
  
   subplot(3,2,2*i-1);
   img = inv_solve(imdl,vv(:,1),vv(:,2:end));
   img.calc_colours.ref_level = 0;
   show_slices(img,[inf,inf,2]);

   subplot(3,2,2*i-0);
   img.elem_data = img.elem_data(:,6);
   if i<3
   show_3d_slices(img,[1.6,2.0,2.4],[],[0.4]); view(-20,20);
   else
   show_slices(img); view(-20,20);
   end
end

print_convert mar29-horse-breathing.jpg

