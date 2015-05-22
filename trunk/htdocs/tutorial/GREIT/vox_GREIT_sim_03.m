%% Set up dual-model
imdl = select_imdl(fmdl,{'Basic GN dif'});
opt = struct();
opt.imgsz = [32 32 3*32];
opt.cube_voxels= true;
opt.prune_model= false; % % GREIT doesn't need a c2f and this is slow
imdl = mk_voxel_volume(imdl,opt);

% prune manually - remove elements with all nodes outside radius
rmdl = imdl.rec_model;
D = sqrt(sum(rmdl.nodes(:,1:2).^2,2)) <= 2;
E = all(D(reshape(rmdl.elems',4*6,[])'),2);
E2= logical(kron(E,true(6,1)));
rmdl.elems = rmdl.elems(E2,:);
rmdl.coarse2fine = rmdl.coarse2fine(E2,E);
rmdl = rmfield(rmdl, 'boundary'); 
rmdl.inside = E2;
imdl.rec_model = rmdl;

%% Set up GREIT
% specify the target distribution
[x, y, z] = ndgrid(linspace(-2,2,20),linspace(-2,2,20),.5:.25:2.5);
idx = ( x(:).^2 + y(:).^2 ) < 1.8^2;
gopt.distr = [x(idx) y(idx) z(idx)]';
% specify noise figure
gopt.noise_figure = 1.0; % if commented out, the given weight is used
% train GREIT
imdl = mk_GREIT_model(imdl,...
                      0.15,... % target radius
                      1,   ... % weight 20.5612 for NF = 1
                      gopt);   % options

%% Reconstruct
rimg = inv_solve(imdl,vh,vi);