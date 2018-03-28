function imdl = solve_RM_2Dslice(imdl, sel_fcn)
% SOLVE_RM_2DSLICE: cut slices out of a 3D model
%  which has reconstruction matrix calculated
%
% Basic usage
% For example, given a 3D GREIT model
%    vopt.zvec= 0.5:0.2:2.5;
%    vopt.imgsz = [32 32];
%    [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
%    imdl = mk_GREIT_model(imdl, 0.20, [], opt);
%
% One could cut the z=1 slice as
%    imdl2= solve_RM_2Dslice(imdl, 1.0);
% 
% Note that the reconstruction planes are between the
% cut planes specified in vopt.zvec
%
% options:
%  imdl.solve_RM_2Dslice.vert_tol = .001;
%       % Vertical tolerance for elems on plane (Default .001)
%  imdl.solve_RM_2Dslice.keep3D = 1; %keep 3D aspect of cuts

% (C) 2015-2018 Andy Adler and Bartlomiej Grychtol
% License: GPL version 2 or version 3
% $Id$

if isstr(imdl) && strcmp(imdl,'UNIT_TEST'); do_unit_test; return; end

% Options
vert_tol = .001;
try
  vert_tol = imdl.solve_RM_2Dslice.vert_tol;
end

ctr = interp_mesh(imdl.rec_model); ct3= ctr(:,3);

if isnumeric( sel_fcn )
   ff = abs(ct3-sel_fcn) < vert_tol;
else
   error huh?
end
imdl.rec_model.coarse2fine(~ff,:) = [];
imdl.rec_model.elems(~ff,:) = [];

fb = any(imdl.rec_model.coarse2fine,1);
imdl.rec_model.coarse2fine(:,~fb) = [];
imdl.fwd_model.coarse2fine(:,~fb) = [];
imdl.solve_use_matrix.RM(~fb,:) = [];
imdl.rec_model.boundary = find_boundary(imdl.rec_model);


imdl.rec_model.nodes(:,3) = [];
imdl.rec_model.elems(:,4) = [];

nelems = imdl.rec_model.elems;
nnodes = unique(nelems(:));
nnidx = zeros(max(nnodes),1);
nnidx(nnodes) = 1:length(nnodes);
nnodes = imdl.rec_model.nodes(nnodes,:);
nelems = reshape(nnidx(nelems),size(nelems));
imdl.rec_model.nodes = nnodes;
imdl.rec_model.elems = nelems;
imdl.rec_model.boundary = find_boundary(imdl.rec_model);

function do_unit_test
   [vh,vi] = test_fwd_solutions;
   % inverse geometry
   fmdl= ng_mk_cyl_models([4,1,.5],[16,1.5,2.5],[0.05]);
   fmdl= mystim_square(fmdl);

   vopt.imgsz = [32 32];
   vopt.zvec = linspace( 0.75,3.25,6);
   vopt.save_memory = 1;
   opt.noise_figure = 2;
   [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
   imdl = mk_GREIT_model(imdl, 0.2, [], opt);

   unit_test_cmp('RM size', size(imdl.solve_use_matrix.RM), [4280,928]);
   unit_test_cmp('RM', imdl.solve_use_matrix.RM(1,1:2), ...
       [7.896475314882005  -3.130412179417599], 1e-10);

   img = inv_solve(imdl, vh, vi);
   unit_test_cmp('img size', size(img.elem_data), [4280,5]);
   [mm,ll] =max(img.elem_data(:,1));
   unit_test_cmp('img', [mm,ll], ...
       [0.426631207850641, 1274], 1e-10);

   img.show_slices.img_cols = 1;
   slice1 = calc_slices(img);
   subplot(131); show_fem(fmdl); title 'fmdl'
   subplot(132); show_slices(img,[inf,inf,2]); title '3D slice';
   imdl = solve_RM_2Dslice(imdl, 2.0);

   unit_test_cmp('RM size', size(imdl.solve_use_matrix.RM), [856,928]);
   unit_test_cmp('RM', imdl.solve_use_matrix.RM(1,1:2), ...
       [-13.546930204647675   9.664897892828284], 1e-10);
   img = inv_solve(imdl, vh, vi);
   unit_test_cmp('img size', size(img.elem_data), [856,5]);
   [mm,ll] =max(img.elem_data(:,1));
   unit_test_cmp('img', [mm,ll], ...
       [0.031608449353163, 453], 1e-10);

   img.show_slices.img_cols = 1;
   subplot(133); show_slices(img); title '2D';
   slice2 = calc_slices(img);

   unit_test_cmp('slice Nan', isnan(slice1), isnan(slice2))
   slice1(isnan(slice1))= 0;
   slice2(isnan(slice2))= 0;
   unit_test_cmp('slice value', slice1, slice2, 1e-10)
   

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
