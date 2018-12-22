function [mapping, outside] = mk_analytic_c2f( f_mdl, c_mdl, opt)
% MK_ANALYTIC_C2F: create a mapping matrix from coarse to fine FEM
% [c2f,out]= mk_analytic_c2f( f_mdl, c_mdl, opt );
%  
% Parameters:
%    c_mdl is coarse fwd_model
%    f_mdl is fine fwd_model
%    opt   is option structure
%
% C2F_ij is the fraction of f_mdl element i which is
%   contained in c_mdl element j. This is used to map
%   from data on the reconstruction model (c_mdl) to
%   the forward model f_mdl as 
%      elem_data_fine = Mapping*elem_data_coase
%
% OUT_i is the fraction of f_mdl element i which is not
%   contained in any c_mdl element.
%
% OPTIONS:
% if the geometry of the fine and coarse models are not
%  aligned, then they can be translated and mapped using
%    coarse_xyz = (M* (fine_xyz - T)')'
%  where
%    T= c_mdl.mk_analytic_c2f.f2c_offset (1xN_dims)
%    M= c_mdl.mk_analytic_c2f.f2c_project (N_dimsxN_dims)
%  by default T= [0,0,0] and M=eye(3)
%
% if c_mdl is 2D and f_mdl is 3D, then parameter
%     c_mdl.mk_analytic_c2f.z_depth
%     indicates the +/- z_depth which elements in 2D are
%     considered to be extruded in 3D (default inf)
%
% See also MK_C2F_CIRC_MAPPING, MK_APPROX_C2F, MK_COARSE_FINE_MAPPING

% (C) 2007-2015 Andy Adler and Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

if ischar(f_mdl) && strcmp(f_mdl, 'UNIT_TEST'); do_unit_test; return; end

if ischar(f_mdl) && strcmp(f_mdl, 'LOAD'); preload; return; end

if nargin== 2, opt = struct; end

opt = assign_defaults( c_mdl, f_mdl, opt );

copt.cache_obj = cache_obj(c_mdl, f_mdl, opt);
copt.fstr = 'mk_analytic_c2f';

[mapping, outside] = eidors_cache(@mapping_calc,{f_mdl,c_mdl,opt},copt);



function [mapping, outside] = mapping_calc(f_mdl, c_mdl, opt)
    
   f_mdl= offset_and_project( f_mdl, opt);
 
   f_dim = elem_dim(f_mdl);
   c_dim = elem_dim(c_mdl);
   d_num = f_dim * 10 + c_dim;
   
   switch d_num
      case 33
         mapping = mk_tet_c2f(f_mdl,c_mdl, opt);
      case 22
         mapping = mk_tri_c2f(f_mdl,c_mdl, opt);
      case 32
         mapping = mk_tri2tet_c2f(f_mdl, c_mdl, opt);
      case 23
         mapping = mk_tri2tet_c2f(c_mdl, f_mdl, opt);
         mapping = bsxfun(@times,   mapping , get_elem_volume(c_mdl));
         if isinf(opt.z_depth)
            height = max(c_mdl.nodes(:,3)) - min(c_mdl.nodes(:,3));
         else
            height = opt.z_depth;
         end
         mapping = mapping * height;
         mapping = bsxfun(@rdivide, mapping', get_elem_volume(f_mdl)); 
      otherwise
         error('EIDORS:WrongDim',['Unsupported combination of dimensions: ' ...
            'c_mdl=%d, f_mdl=%d'],c_dim,f_dim);
   end

   
   if isfield(c_mdl,'coarse2fine')
      mapping = mapping*c_mdl.coarse2fine;
   end
   
   if nargout>1;
      outside = 1 - sum(mapping,2);
   end


% Mapping depends only on nodes and elems - remove the other stuff
function c_obj = cache_obj(c_mdl, f_mdl, opt)
   c_obj = {c_mdl.nodes, c_mdl.elems,  ...
            f_mdl.nodes, f_mdl.elems, opt};

% Offset and project f_mdl as required
function f_mdl= offset_and_project( f_mdl, opt)
    [fn,fd]= size(f_mdl.nodes);
    T= opt.f2c_offset;
    M= opt.f2c_project;
    
    f_mdl.nodes= (M*( f_mdl.nodes - ones(fn,1)*T )')';

function opt = assign_defaults( c_mdl, f_mdl, org )
   opt = struct;
   try opt = c_mdl.mk_analytic_c2f; end
   fn = fieldnames(org);
   for f = 1:numel(fn)
      opt.(fn{f}) = org.(fn{f}); % explicit options override those on the model
   end
   [fn,fd]= size(f_mdl.nodes);
   try    opt.f2c_offset; % test exist
   catch  opt.f2c_offset= zeros(1,fd);
   end
   try    opt.f2c_project;
   catch  opt.f2c_project= speye(fd);
   end
   try    opt.z_depth;
   catch  opt.z_depth= inf;
   end

     
    
function do_unit_test
   ll = eidors_msg('log_level', 1);
   do_original_2d_tests
   test_2d3d_handling
   mk_tri_c2f UNIT_TEST
   mk_tet_c2f UNIT_TEST
   mk_tri2tet_c2f UNIT_TEST
   eidors_msg('log_level',ll);
   
function test_2d3d_handling
   tet = mk_grid_model([],[0 1],[0 1],[0 1]);
   tet = rmfield(tet,'coarse2fine');
   tri = mk_grid_model([],[0 1],[0 1]);
   tri = rmfield(tri,'coarse2fine');
   
   c2f_a = mk_analytic_c2f(tri,tet);
   c2f_a = bsxfun(@times, c2f_a, get_elem_volume(tri));
   c2f_b = mk_analytic_c2f(tet,tri);
   c2f_b = bsxfun(@times, c2f_b, get_elem_volume(tet));
   
   unit_test_cmp('2d3d:01',sum(c2f_a(:)), 1,eps);
   unit_test_cmp('2d3d:02',sum(c2f_b(:)), 1,eps);
   unit_test_cmp('2d3d:03',c2f_a,c2f_b',eps)
   

   
function do_original_2d_tests
   fmdl = mk_circ_tank(2,[],2); fmdl.nodes = fmdl.nodes*2;
   cmdl = mk_circ_tank(2,[],2); cmdl.nodes = cmdl.nodes*2;
   c2f = mk_analytic_c2f( fmdl, cmdl);
   unit_test_cmp('t1',c2f,eye(16),eps)

   fmdl = mk_circ_tank(3,[],2);
   fmdl.nodes = fmdl.nodes*3;
   c2f = mk_analytic_c2f( fmdl, cmdl);
   unit_test_cmp('t2',c2f,[eye(16);zeros(20,16)],eps)

   fmdl = mk_circ_tank(2,[],2); fmdl.nodes = fmdl.nodes*2;
   cmdl = mk_circ_tank(1,[],2); cmdl.nodes = cmdl.nodes*1;
   c2f = mk_analytic_c2f( fmdl, cmdl);
   unit_test_cmp('t3',c2f,[eye(4);zeros(12,4)],eps)

   cmdl = mk_circ_tank(1,[],2); cmdl.nodes = cmdl.nodes*0.8;
   c2f = mk_analytic_c2f( fmdl, cmdl);
   unit_test_cmp('t4',c2f,[eye(4)*.64;zeros(12,4)], eps)


   fmdl = mk_circ_tank(10,[],2);
   cmdl = mk_circ_tank(8,[],2);
   c2f = mk_analytic_c2f( fmdl, cmdl);
   unit_test_cmp('t6',sum(c2f(1:324,:)'),ones(1,324),1e-14);

   cmdl.nodes = cmdl.nodes*0.95;
   % show_fem(fmdl); hold on ; show_fem(cmdl); hold off
   c2f = mk_analytic_c2f( fmdl, cmdl);

% preload into cache before tests
function preload

   % Create forward, fine tank model
   electrodes_per_plane = 16;
   number_of_planes = 2;
   tank_radius = 0.2;
   tank_height = 0.5;
   fine_mdl = ng_mk_cyl_models([tank_height,tank_radius],...
       [electrodes_per_plane,0.15,0.35],[0.01]);
    
   % Create coarse model for inverse problem
   coarse_mdl_maxh = 0.07; % maximum element size 
   coarse_mdl = ng_mk_cyl_models([tank_height,tank_radius,coarse_mdl_maxh],[0],[]);

   disp('Calculating coarse2fine mapping ...');
   inv3d.fwd_model.coarse2fine = ...
          mk_analytic_c2f( fine_mdl, coarse_mdl);
   disp('   ... done');
