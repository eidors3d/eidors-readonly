function [mapping, outside] = mk_coarse_fine_mapping( f_mdl, c_mdl );
% MK_COARSE_FINE_MAPPING: create a mapping matrix from coarse to fine FEM
% [c2f,out]= mk_coarse_fine_mapping( f_mdl, c_mdl );
%  
% Parameters:
%    c_mdl is coarse fwd_model
%    f_mdl is fine fwd_model
%
% C2F_ij is the fraction if f_mdl element i which is
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
%    coarse_xyz = M*( fine_xyz - T)
%  where
%    T= c_mdl.mk_coarse_fine_mapping.f2c_offset (1xN_dims)
%    M= c_mdl.mk_coarse_fine_mapping.f2c_project (N_dimsxN_dims)
%  by default T= [0,0,0] and M=eye(3)
%
% if c_mdl is 2D and f_mdl is 3D, then parameter
%     c_mdl.mk_coarse_fine_mapping.z_depth
%     indicates the +/- z_depth which elements in 2D are
%     considered to be extruded in 3D (default inf)
%
% NOTES:
% if c_mdl and f_mdl do not cover the same area, then 
%    sum(c2f') will not be 1. If all coarse elements cover
%    at least partially the fine ones, then this 
%    can be corrected by:
%      c2f = c2f./(sum(c2f,2) * ones(1,size(c2f,2)));
% if not all coarse elements cover fine ones, then this
%    approach cannot be used. This will be fixed in a 
%    future release

%
% See also MK_C2F_CIRC_MAPPING

% (C) 2007-2012 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(f_mdl) && strcmp(f_mdl, 'UNIT_TEST'); do_unit_test; return; end

if isstr(f_mdl) && strcmp(f_mdl, 'LOAD'); load; return; end

[c_mdl, f_mdl] = assign_defaults( c_mdl, f_mdl );
c_obj = cache_obj(c_mdl, f_mdl);

f_mdl= offset_and_project( f_mdl, c_mdl);
mapping = eidors_obj('get-cache', c_obj, 'coarse_fine_mapping');
if ~isempty(mapping)
    eidors_msg('mk_coarse_fine_mapping: using cached value', 3);
else

    z_depth = c_mdl.mk_coarse_fine_mapping.z_depth;

    f_elems = all_contained_elems( f_mdl, c_mdl, z_depth);
    mapping = contained_elems_i( f_mdl, c_mdl, f_elems, z_depth);

    if isfield(c_mdl,'coarse2fine')
       mapping = mapping*c_mdl.coarse2fine;
    end

    eidors_obj('set-cache', c_obj, 'coarse_fine_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end

if nargout>1;
  outside = 1 - sum(mapping,2);
end

% Mapping depends only on nodes and elems - remove the other stuff
function c_obj = cache_obj(c_mdl, f_mdl)
   c_obj = {c_mdl.nodes, c_mdl.elems, c_mdl.mk_coarse_fine_mapping, ...
            f_mdl.nodes, f_mdl.elems, f_mdl.interp_mesh};


% find all elems of ff_mdl completely contained in cc_mdl
function c_elems = all_contained_elems( fm, cm, z_depth)
    [nf,ef]= size(fm.elems);
    [nc,ec]= size(cm.elems);
    fm_pts = zeros(nf*ef,3);
    % shrink pts slightly (s_fac) so they're not on the boundary
    % by shrinking, we avoid cases where an element is
    % only slighly intersecting another. This is beyond the
    % resolution of the next step (interpolation) anyway
    s_fac= .9; % .9999;
    for dim= 1:ef-1 % ef-1 is dimensions in fm
       % fm_pts is local_nodes x elems x xyz
       fm_pt= reshape(fm.nodes(fm.elems,dim),nf,ef);
       fm_ctr= mean(fm_pt,2)*ones(1,ef);
       fm_pt = s_fac*(fm_pt-fm_ctr) + fm_ctr;
       fm_pts(:,dim) = fm_pt(:);
    end

    tsn= search_fm_pts_in_cm(cm, fm_pts, z_depth);
    tsn= reshape( tsn, [], ef);
    % if all points are outside (NaN) then c_elems = -1
    % if all points are in one elem   then c_elems = elem #
    % if all points are in diff elems then c_elems = 0
    c_elems= all(diff(tsn,1,2)==0,2) .* tsn(:,1);
    c_elems(all(tsn==-1,2))= -1; % all points outside


% tsn = vector of length size(fm_pts,1)
% tsn(i) = elem in cm which contains point
% tsn(i) = -1 if point is outside cm (and z_depth, if appropriate)
function tsn= search_fm_pts_in_cm(cm, fm_pts, z_depth);
    dc= size(cm.elems,2)-1;  %coarse dim
    df= size(fm_pts,2); %fine dim

    tsn= -ones(size(fm_pts,1),1);
    not_oor= (tsn==-1); % logical 1

    if dc==2  %corse model is 2D

       if df==3
       % look for f_mdl z not out of range 
          not_oor= not_oor &  any( abs(fm_pts(:,3) ) <= z_depth , 2);
       end
       dims=1:2;

    elseif dc==3 %coarse model is 3D

       dims=1:3; 

    else
       error('coarse model must be 2 or 3D');
    end

    tsn(not_oor)= tsearchn(cm.nodes(:,dims), cm.elems, fm_pts(not_oor,dims));
    tsn(isnan(tsn))= -1;

% find fraction of elem contained in cm
% idx is the index into which the elem is contained
% if idx >= 1, the element is completely with in that coarse elem
% if idx == 0, the element is crosses several coarse elems
% if idx ==-1, the element is outside the coarse model
function c_elems = contained_elems_i( fm, cm, idx, z_depth)
   [nc,dc]= size(cm.elems);
   [nf,df]= size(fm.elems);

   fidx= find(idx==0);
   l_fidx= length(fidx);

   c_e_i= []; c_e_j=[]; c_e_v=[];

   % lower density interpolation in higher dimentions, since 
   % the added dimensions will give extra interpolation points.
   n_interp = 7-df;
   interp_mdl.nodes= fm.nodes;
   interp_mdl.elems= fm.elems(fidx,:);
   interp_mdl.interp_mesh.n_interp = fm.interp_mesh.n_interp;
 

   fm_pts = interp_mesh( interp_mdl, n_interp);
   l_interp = size(fm_pts,3);

   fm_pts = permute( fm_pts, [3,1,2]);
   fm_pts = reshape(fm_pts, [], df-1);

   tsn= search_fm_pts_in_cm(cm, fm_pts, z_depth);

   tsn_idx= ones(l_interp,1)*fidx(:)';
   tsn_idx= tsn_idx(:);
   % find and isolate outside elements
   outside_idx= tsn==-1;
   tsn(outside_idx) = [];
   tsn_idx(outside_idx) = [];
   
   in_idx= find(idx<=0);
   ridx= 1:nf; ridx(in_idx)=[];
   idx(in_idx)=[];

   % first term is contribution from f_elems in one c_elem
   % next term is contribution from f_elems in many c_elems, weighted
   c_elems = sparse(ridx,idx,1,nf,nc) +  ...
             sparse(tsn_idx,tsn,1,nf,nc)/l_interp;
      
% Do 3D interpolation of region xyzmin= [x,y,z] to xyzmax
%  with n_interp points in the minimum direction
function xyz = interpxyz( xyzmin, xyzmax, n_interp)
    xyzdelta= xyzmax - xyzmin;
    xyz_interp = 1 + floor(n_interp * xyzdelta / max(xyzdelta) );
    xspace = linspace(xyzmin(1), xyzmax(1), xyz_interp(1) );
    yspace = linspace(xyzmin(2), xyzmax(2), xyz_interp(2) );
    zspace = linspace(xyzmin(3), xyzmax(3), xyz_interp(3) );
    [xx3,yy3,zz3] = ndgrid( xspace, yspace, zspace );
    xyz= [xx3(:), yy3(:), zz3(:)];

% Offset and project f_mdl as required
function f_mdl= offset_and_project( f_mdl, c_mdl)
    [fn,fd]= size(f_mdl.nodes);
    T= c_mdl.mk_coarse_fine_mapping.f2c_offset;
    M= c_mdl.mk_coarse_fine_mapping.f2c_project;

    f_mdl.nodes= ( f_mdl.nodes - ones(fn,1)*T )*M;

function [c_mdl f_mdl] = assign_defaults( c_mdl, f_mdl )
    [fn,fd]= size(f_mdl.nodes);
    try    c_mdl.mk_coarse_fine_mapping.f2c_offset; % test exist
    catch  c_mdl.mk_coarse_fine_mapping.f2c_offset= zeros(1,fd);
    end
    try    c_mdl.mk_coarse_fine_mapping.f2c_project;
    catch  c_mdl.mk_coarse_fine_mapping.f2c_project= speye(fd);
    end
    try    c_mdl.mk_coarse_fine_mapping.z_depth;
    catch  c_mdl.mk_coarse_fine_mapping.z_depth= inf;
    end
    try    f_mdl.interp_mesh.n_interp;
    % lower density interpolation in higher dimentions, since
    % the added dimensions will give extra interpolation points.
    catch  f_mdl.interp_mesh.n_interp = 7 - size(f_mdl.elems,2);
    end
     
    
function do_unit_test
    fmdl = mk_circ_tank(2,[],2); fmdl.nodes = fmdl.nodes*2;
    cmdl = mk_circ_tank(2,[],2); cmdl.nodes = cmdl.nodes*2;
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);
    unit_test_cmp('t1',c2f,eye(16))

    fmdl = mk_circ_tank(3,[],2);
    fmdl.nodes = fmdl.nodes*3;
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);
    unit_test_cmp('t2',c2f,[eye(16);zeros(20,16)])

    fmdl = mk_circ_tank(2,[],2); fmdl.nodes = fmdl.nodes*2;
    cmdl = mk_circ_tank(1,[],2); cmdl.nodes = cmdl.nodes*1;
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);
    unit_test_cmp('t3',c2f,[eye(4);zeros(12,4)])

    cmdl = mk_circ_tank(1,[],2); cmdl.nodes = cmdl.nodes*0.8;
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);
    unit_test_cmp('t3',c2f,[eye(4)*2/3;zeros(12,4)])

    cmdl = mk_circ_tank(1,[],2); cmdl.nodes = cmdl.nodes*1.2;
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);
    unit_test_cmp('t3',c2f,[eye(4);eye(4)/3;kron(eye(4),[1;1])/15]);

    fmdl = mk_circ_tank(10,[],2);
    cmdl = mk_circ_tank(8,[],2);
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);
    unit_test_cmp('t4',sum(c2f'),ones(1,size(c2f,1)),1e-14);
   
    cmdl.nodes = cmdl.nodes*0.95;
% show_fem(fmdl); hold on ; show_fem(cmdl); hold off
    c2f = mk_coarse_fine_mapping( fmdl, cmdl);

function load

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
       mk_coarse_fine_mapping( fine_mdl, coarse_mdl);
disp('   ... done');
