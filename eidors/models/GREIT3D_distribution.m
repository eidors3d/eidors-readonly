function [imdl,distr] = GREIT3D_distribution(fmdl, vopt) 
% GREIT3D_distribution: create target distributions for 3D GREIT
%
% The basic usage is part of creating a 3D GREIT model, 
%  given a forward model fmdl
%    vopt.imgsz = [32 32 19]; % and other parameters
%    opt.noise_figure = NF;   % and other GREIT options
%    [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
%    imdl = mk_GREIT_model(imdl, 0.20, [], opt);
%
% INPUTS:
%  fmdl => EIDORS forward model structure
%  vopt => Vertical options (see mk_voxel_volume help)
%
% Options for vopt
%   1: Simple
%        vopt.imgsz = [32 32 19];
%        vopt.square_pixels = true;
%   2: Specify in detail the z planes (e.g.)
%        vopt.imgsz = [32 32];
%        vopt.zvec = linspace( 0.5,3.5,15);
%          (Vertical size is: length(vopt.zvec)-1
%   3: Specify x,y and z planes (e.g.)
%        vopt.xvec = linspace(-.95,.95,32);
%        vopt.yvec = linspace(-.95,.95,32);
%        vopt.zvec = linspace( .1,3.9,11);
%
% NOTE: mk_voxel_volume is slow and uses lots of memory
%   To reduce the memory footprint, try
%   vopt.save_memory = 1; %(or try 10)

% (C) 2015-2018 Andy Adler and Bartlomiej Grychtol
% License: GPL version 2 or version 3
% $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

   imdl = select_imdl(fmdl,{'Basic GN dif'});
   imdl = mk_voxel_volume(imdl,vopt);
   msm = imdl.rec_model.mdl_slice_mapper;
   imdl.rec_model.mdl_slice_mapper.y_pts = ...
       fliplr(imdl.rec_model.mdl_slice_mapper.y_pts);
   [x, y, z] = ndgrid(msm.x_pts, msm.y_pts, msm.z_pts);

   IN = point_in_vox(imdl.rec_model,[x(:),y(:),z(:)]);
   distr = [x(IN),y(IN),z(IN)]';

function out = point_in_vox(fmdl,points, epsilon)

   if nargin < 3
      epsilon = 0;
   end

   levels = unique(fmdl.nodes(:,3));
   out = false(numel(points)/3,1);

   opt.faces = true;
   fmdl = fix_model(fmdl, opt);

   jnk.nodes = fmdl.nodes(:,1:2);
   jnk.type = 'fwd_model';

   for i = 1:numel(levels)-1
      idx = (points(:,3) - levels(i)) > epsilon;
      idx = idx & ((levels(i+1) - points(:,3)) > epsilon);
      
      nidx = fmdl.nodes(:,3) == levels(i);
      F = all(nidx(fmdl.faces),2);
      jnk.elems = fmdl.faces(F,:);
      bnd = find_boundary(jnk);
      bnd = order_loop(jnk.nodes(unique(bnd(:)),:));
      out(idx) = out(idx) | inpolygon(points(idx,1),points(idx,2),bnd(:,1),bnd(:,2));
   end

function do_unit_test
   [vh,vi] = test_fwd_solutions;
   % inverse geometry
   fmdl= ng_mk_cyl_models([4,1,.5],[16,1.5,2.5],[0.05]);
   fmdl= mystim_square(fmdl);

   vopt.imgsz = [32 32];
   vopt.zvec = linspace( 0.5,3.5,7);
   vopt.save_memory = 1;
%  opt.noise_figure = 2;
   [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
   n_weight = 18.1608; % approx NF=2
   imdl = mk_GREIT_model(imdl, 0.2, n_weight, opt);

   unit_test_cmp('RM size', size(imdl.solve_use_matrix.RM), [5136,928]);
   tst = imdl.solve_use_matrix.RM(1,1:2);
   tst_= [8.572908305490031 -18.600780994066596];
   unit_test_cmp('RM', tst, tst_, 1e-10);

   img = inv_solve(imdl, vh, vi);
   unit_test_cmp('img size', size(img.elem_data), [5136,5]);
   [mm,ll] =max(img.elem_data(:,1));
   unit_test_cmp('img', [mm,ll], ...
       [0.582158646727019, 1308], 1e-10);

   img.show_slices.img_cols = 1;
   subplot(131); show_fem(fmdl); title 'fmdl'
   subplot(132); show_slices(img,[inf,inf,2]); title 'centre';
   subplot(133); show_slices(img,[inf,inf,1.5]); title 'bottom';
   

function fmdl = mystim_square(fmdl);
   [~,fmdl] = elec_rearrange([16,2],'square', fmdl);
   [fmdl.stimulation, fmdl.meas_select] =  ...
       mk_stim_patterns(32,1,[0,5],[0,5],{},1);

function [vh,vi] = test_fwd_solutions;
   posns= linspace(1.0,3.0,5);
   str=''; for i=1:length(posns);
      extra{i} = sprintf('ball%03d',round(posns(i)*100));
      str = [str,sprintf('solid %s=sphere(0.5,0,%f;0.1); ', extra{i}, posns(i))];
   end
   extra{i+1} = str;
   fmdl= ng_mk_cyl_models([4,1,.2],[16,1.5,2.5],[0.05],extra); 
   fmdl = mystim_square(fmdl);
   
   img= mk_image(fmdl,1);
   vh = fwd_solve(img);
   %vh = add_noise(2000, vh);
   for i=1:length(posns);
      img= mk_image(fmdl,1);
      img.elem_data(fmdl.mat_idx{i+1}) = 2;
      vi{i} = fwd_solve(img);
   end;
   vi = [vi{:}];
