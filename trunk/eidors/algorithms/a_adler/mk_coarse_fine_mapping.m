function mapping = mk_coarse_fine_mapping( f_mdl, c_mdl );
% MK_COARSE_FINE_MAPPING: create a mapping matrix from coarse to fine FEM
% mapping= mk_coarse_fine_mapping( f_mdl, c_mdl );
%
% Mapping approximates elem_data_fine from elem_data_coase
%   elem_data_fine = Mapping*elem_data_coase
%
% c_mdl is coarse fwd_model
% f_mdl is fine fwd_model
%
% if the geometry of the fine and coarse models are not
%  aligned, then they can be translated and mapped using
%    coarse_xyz = M*( fine_xyz - T)
%  where
%    T= c_mdl.mk_coarse_fine_mapping.f2c_offset (1xN_dims)
%    M= c_mdl.mk_coarse_fine_mapping.f2c_project (N_dimsxN_dims)
%  by default T= [0,0,0] and M=eye(3)
%
% if c_mdl is 2D and f_mdl is 3D, then parameter
%     c_mdl.mk_coarse_fine_mapping.z_depth (default = inf)
%     indicates the +/- z_depth which elements in 2D are
%     considered to be extruded in 3D

% (C) 2007-2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

c_obj = cache_obj(c_mdl, f_mdl);

f_mdl= offset_and_project( f_mdl, c_mdl);
mapping = eidors_obj('get-cache', c_obj, 'coarse_fine_mapping');
if ~isempty(mapping)
    eidors_msg('mk_coarse_fine_mapping: using cached value', 3);
else

    try
       z_depth = c_mdl.mk_coarse_fine_mapping.z_depth;
    catch
       z_depth = inf;
    end

    f_elems = all_contained_elems( f_mdl, c_mdl, z_depth);
    mapping = contained_elems_i( f_mdl, c_mdl, f_elems, z_depth);

    if isfield(c_mdl,'coarse2fine')
       mapping = mapping*c_mdl.coarse2fine;
    end

    eidors_obj('set-cache', c_obj, 'coarse_fine_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end

% Mapping depends only on nodes and elems - remove the other stuff
function c_obj = cache_obj(c_mdl, f_mdl)
   c_obj = {c_mdl.nodes, c_mdl.elems,  ...
            f_mdl.nodes, f_mdl.elems};

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
function xyz = interpxyz( xyzmin, xyzmax, n_interp);
    xyzdelta= xyzmax - xyzmin;
    xyz_interp = 1 + floor(n_interp * xyzdelta / max(xyzdelta) );
    xspace = linspace(xyzmin(1), xyzmax(1), xyz_interp(1) );
    yspace = linspace(xyzmin(2), xyzmax(2), xyz_interp(2) );
    zspace = linspace(xyzmin(3), xyzmax(3), xyz_interp(3) );
    [xx3,yy3,zz3] = ndgrid( xspace, yspace, zspace );
    xyz= [xx3(:), yy3(:), zz3(:)];

% Offset and project f_mdl as required
function f_mdl= offset_and_project( f_mdl, c_mdl);
    [fn,fd]= size(f_mdl.nodes);
    try
       T= c_mdl.mk_coarse_fine_mapping.f2c_offset;
    catch
       T= zeros(1,fd);
    end
    try
       M= c_mdl.mk_coarse_fine_mapping.f2c_project;
    catch
       M= speye(fd);
    end

    f_mdl.nodes= ( f_mdl.nodes - ones(fn,1)*T )*M;
