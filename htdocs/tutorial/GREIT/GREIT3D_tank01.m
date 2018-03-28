%STEP 1: Simple model
[stim,msel] = mk_stim_patterns(32,1,[0,5],[0,5],{},1);
posns= linspace(0.5,3.5,13);
str=''; for i=1:length(posns);
   extra{i} = sprintf('ball%03d',round(posns(i)*100));
   str = [str,sprintf('solid %s=sphere(0.5,0,%f;0.1); ', extra{i}, posns(i))];
end
extra{i+1} = str;
fmdl= ng_mk_cyl_models([4,1,.2],[16,1.5,2.5],[0.05],extra); 
fmdl = mdl_normalize(fmdl, 0);
[~,fmdl] = elec_rearrange([16,2],'square', fmdl);
fmdl.stimulation= stim; fmdl.meas_select = msel;

img= mk_image(fmdl,1);
img.elem_data(vertcat(fmdl.mat_idx{2:end}))= 2;
opt.viewpoint = struct('az',0,'el',10);
show_fem_enhanced(img,opt);

% STEP 2: Simulate data
img= mk_image(fmdl,1);
vh = fwd_solve(img);
for i=1:length(posns)-4;
   img= mk_image(fmdl,1);
   img.elem_data(fmdl.mat_idx{i+1}) = 2;
   vi{i} = fwd_solve(img);
end;

clear opt;

% STEP 3: Reconstruction model
fmdl= ng_mk_cyl_models([4,1,.5],[16,1.5,2.5],[0.05]);
fmdl = mdl_normalize(fmdl, 0);
[~,fmdl] = elec_rearrange([16,2],'square', fmdl);
fmdl.stimulation= stim; fmdl.meas_select = msel;

% STEP 4: Reconstruction using GREIT 2D
   opt.imgsz = [32 32];
   opt.square_pixels = true;
   opt.noise_figure = 0.5;
   img = mk_image(fmdl,1);
   imdl2 = mk_GREIT_model(img, 0.25, [], opt);

subplot(311);
img = inv_solve(imdl2,vh,[vi{:}]);
img.show_slices.img_cols= 9;
show_slices(img);

% STEP 4: Reconstruction using GREIT 3D
   vopt.imgsz = [32 32];
   vopt.zvec = linspace( 0.75,3.25,6);
   vopt.save_memory = 1;
   opt.noise_figure = 2.0;
   [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
   imdl3= mk_GREIT_model(imdl, 0.20, [], opt);
subplot(312);
img = inv_solve(imdl3,vh,[vi{:}]);
img.show_slices.img_cols= 9;
show_slices(img,[inf,inf,1.5;
                 inf,inf,2;]);

subplot(313);
imdl3a= solve_RM_2Dslice(imdl3,1.5);
img = inv_solve(imdl3a,vh,[vi{:}]);
img.show_slices.img_cols= 9;
show_slices(img);
