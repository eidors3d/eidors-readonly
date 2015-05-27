function stim_patterns_sim
for i=1:3; switch i
      case 1; sel = 0; weight = 2; ofs = 0;
      case 2; sel = 1; weight = 2; ofs = 3;
      case 3; sel = 2; weight = 2; ofs = 6;
   end;

   fmdl = mk_forward_model( sel )
   [imgi,vh,vi] = mk_sims(fmdl);
%  show_image(imgi, 0);

   imdl = GREIT_3D_inv(fmdl, weight);
   imgr = inv_solve(imdl,vh,vi);
   imgr.calc_colours.ref_level = 0;
   show_image(imgr, ofs);
end


opt.horz_cut     = 20;
opt.vert_cut     = 20;
print_convert('greit-3d-example.jpg',opt);

function fmdl = mk_forward_model( sel )
%% Forward model
   fmdl= ng_mk_cyl_models([3,2,.4],[16,1,2],[.1,0,.025]);
   fmdl.stimulation = mk_stim_patterns(32,1,[0,5],[0,5],{},1);

   idx = 0:31; i1 = floor(idx/16); i2= rem(idx,16);
   switch sel
      case 0; % "Planar"
      case 1; % "Zigzag offset"
          idx = i2*2 + i1;
      case 2; % "square"
          idx = i2*2 + i1.*(rem(i2,2)==0) + (~i1).*(rem(i2,2)==1);;
      otherwise; error('huh?');
   end

%  show= 1:32; show=show(idx+1); disp(reshape(show,16,2)');
   fmdl.electrode(idx+1) = fmdl.electrode;
   %show_fem(fmdl,[0,1]);

function [imgi,vh,vi] = mk_sims(fmdl);
   % homogeneous
   imgh = mk_image(fmdl,1);
   vh = fwd_solve(imgh);
   % spherical target
   select_fun = inline([ ... 
        '((x-.6).^2+(y-.5).^2+(z-1.05).^2<=0.2^2) |', ...
        '((x+.6).^2+(y-.5).^2+(z-1.95).^2<=0.2^2)'], 'x','y','z');
   % inhomogeneous
   imgi = imgh;
   imgi.elem_data = imgh.elem_data - 0.1*elem_select(fmdl,select_fun);
   vi = fwd_solve(imgi);
   rng('default');
   vi = add_noise(20,vi,vh);


function show_image(img, ofs);
   subplot(3,3,1+ofs)
   show_fem(img,[0,1]);
    view(0,22)

   subplot(3,3,2+ofs)
   show_3d_slices(img,1.25,.5,.5);% cuts through the target center
    view(0,22)

   subplot(3,3,3+ofs)
   levels(:,3) = 2.25:-.5:.25; % cuts through the target center
   levels(:,1:2) = Inf;
   show_slices(img,levels)

function imdl = GREIT_3D_inv( fmdl, weight);
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
   %gopt.noise_figure = 2.0; % if commented out, the given weight is used
   % train GREIT
   imdl = mk_GREIT_model(imdl,...
                         0.15,... % target radius
                         weight,   ... % weight 20.5612 for NF = 1
                         gopt);   % options


