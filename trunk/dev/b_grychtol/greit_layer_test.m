n_elecs = 8;
layers = [2.8,3.2];
stim =  mk_stim_patterns(n_elecs*length(layers),1,[0,1],[0,1], ...
             {'no_meas_current'}, 1);

extra={'lungs','solid lungs = sphere(0.9,0.1,3.65;0.3) or sphere(-0.9,0.1,2.35;0.3);'};
[fmdl,midx] = ng_mk_cyl_models([6, 2,0.1] ,[n_elecs,layers],[0.1 0 0.01], extra);
fmdl.stimulation =  stim;

img = mk_image(fmdl,1); % Homogeneous background
vh = fwd_solve(img);
img.elem_data(midx{2}) = 0.5; % Lung regions
vi = fwd_solve(img);
subplot(121)
show_fem(img,[0,1]); view(0,70);
view(0,10);

opt.imgsz = [32 32];
opt.distr = 3; % non-random, uniform
opt.Nsim = 1000;
opt.target_size = 0.05; % Target size (frac of medium)
% opt.noise_figure = 1.0; % Recommended NF=0.5;

%%
fmdl = ng_mk_cyl_models([6, 2] ,[n_elecs,layers],[0.1]);
fmdl.stimulation =  stim;
fmdl = mdl_normalize(fmdl, 0);
img = mk_image(fmdl,1);


opt.target_plane = [2.5 3 3.5];

imdl = mk_GREIT_model(img, 0.25, 0.3, opt);
rimg = inv_solve(imdl,vh,vi);
subplot(122)
show_fem(fmdl);
hold on
rimg.fwd_model.boundary = find_boundary(rimg.fwd_model);
hh = show_fem(rimg);
set(hh,'linewidth',2)

hold off
view(0,10);
return
%%
tl = linspace(0.5,1.5,3);
opt.target_size = 0.005; % Target size (frac of medium)
for i=1:length(tl);
   opt.target_plane =  tl(i);
   opt.z_depth = 0.2;
   imdl = mk_GREIT_model(img, 0.25, [], opt);
   imdl.target_plane = tl(i);
   imdls(i) = imdl;
end


for i=1:length(tl);
   rimg = inv_solve(imdls(i),vh,vi); %Reconstruct
   rimg.calc_colours.ref_level=0;
   rimgs(i) = rimg;
end

show_slices(rimgs); axis square;