%% Forward model
fmdl= ng_mk_cyl_models([3,2,.4],[16,1,2],[.1,0,.025]);
show_fem(fmdl);
fmdl.stimulation = mk_stim_patterns(16,2,'{ad}','{ad}');

%% Simulation
imgh = mk_image(fmdl,1);
vh = fwd_solve(imgh);
imgi = imgh;
select_fun = inline('(x-.5).^2+(y-.5).^2+(z-1.25).^2<=0.1^2','x','y','z');
imgi.elem_data = imgh.elem_data + elem_select(fmdl,select_fun);
vi = fwd_solve(imgi);

%% Reconstruction
% build a GN inverse model
imdl = select_imdl(fmdl,{'Basic GN dif'});
% prepare dual model
opt.imgsz = [32 32];
opt.zvec = [0 .75:.5:2.25 3];
imdl = mk_voxel_volume(imdl,opt);
% reconstruct
imdl.hyperparameter.value = 1e-7;
rimg = inv_solve(imdl,vh,vi);
show_fem(rimg);
