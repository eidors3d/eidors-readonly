fmdl= ng_mk_cyl_models([3,2,.4],[16,1,2],[.1,0,.025]);
show_fem(fmdl);
fmdl.stimulation = mk_stim_patterns(16,2,'{ad}','{ad}');
imdl = select_imdl(fmdl,{'Basic GN dif'});
%%
opt.imgsz = [32 32 3*32];
% opt.zvec = [0 .75:.5:2.25 3];
opt.cube_voxels= true;
opt.prune_model= false;
imdl = mk_voxel_volume(imdl,opt);

% prune manually
rmdl = imdl.rec_model;
D = sqrt(sum(rmdl.nodes(:,1:2).^2,2)) <= 2;
E = all(D(reshape(rmdl.elems',4*6,[])'),2);
E2= logical(kron(E,true(6,1)));
rmdl.elems = rmdl.elems(E2,:);
rmdl.coarse2fine = rmdl.coarse2fine(E2,E);
rmdl = rmfield(rmdl, 'boundary');
rmdl.inside = E2;
imdl.rec_model = rmdl;

%%
[x, y, z] = ndgrid(linspace(-2,2,20),linspace(-2,2,20),.5:.25:2.5);
idx = ( x(:).^2 + y(:).^2 ) < 1.8^2;
gopt.distr = [x(idx) y(idx) z(idx)]';
gopt.desired_solution_fn = @new_GREIT_desired_soln;
gopt.noise_figure = 0.5;
% gopt.target_size
imdl = mk_GREIT_model(imdl,0.15,0.2, gopt);
%%
imgh = mk_image(fmdl,1);
vh = fwd_solve(imgh);
%%
imgi = imgh;
select_fun = inline('(x-.5).^2+(y-.5).^2+(z-1).^2<=0.1^2','x','y','z');
imgi.elem_data = imgh.elem_data + elem_select(fmdl,select_fun);
vi = fwd_solve(imgi);
%%
rimg = inv_solve(imdl,vh,vi);
show_fem(rimg);

% TODO: desired solution cannot assume equally spaced points