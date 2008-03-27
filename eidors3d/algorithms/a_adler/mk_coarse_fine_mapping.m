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
% $Id: mk_coarse_fine_mapping.m,v 1.20 2008-03-27 17:15:47 aadler Exp $

% Mapping depends only on nodes and elems - remove the other stuff
try; c_mdl= rmfield(c_mdl,'electrode');   end
try; c_mdl= rmfield(c_mdl,'stimulation'); end
try; f_mdl= rmfield(f_mdl,'electrode');   end
try; f_mdl= rmfield(f_mdl,'stimulation'); end

mapping = eidors_obj('get-cache', {f_mdl,c_mdl}, 'coarse_fine_mapping');
if ~isempty(mapping)
    eidors_msg('mk_coarse_fine_mapping: using cached value', 3);
else
    f_mdl= offset_and_project( f_mdl, c_mdl);

    try
       z_depth = c_mdl.mk_coarse_fine_mapping.z_depth;
    catch
       z_depth = inf;
    end

    f_elems = all_contained_elems( f_mdl, c_mdl, z_depth);
    mapping = contained_elems_i( f_mdl, c_mdl, f_elems);

    eidors_obj('set-cache', {f_mdl,c_mdl}, 'coarse_fine_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end


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
    % if all points are outside (NaN) then c_elems = NaN
    % if all points are in one elem   then c_elems = 1
    % if all points are in diff elems then c_elems = 0
    c_elems= all(diff(tsn,1,2)==0,2);
    c_elems(any(isnan(tsn),1))= NaN;

return
    % This is not quite correct, 

    % if node in fm is outside cm, then set it back inside
    % this way only elems in fm that are completely outside cm
    % will be identified that way
    tsn= sort( tsn, 2); % send Nans to right
    isn= find(isnan(tsn));
    tsn(isn)= tsn( 1+rem(isn-1,nf ) );
    c_elems= c_elems.* tsn(:,1);

function tsn= search_fm_pts_in_cm(cm, fm_pts, z_depth);
    dc= size(cm.elems,2);  %coarse dim+1
    [nf,df]= size(fm_pts); %fine dim+1

    tsn= NaN*ones(size(fm_pts,1),1);
    not_oor= ~any( isnan(fm_pts,2) );

    if dc==3
       if df==4
       % look for f_mdl z not out of range 
          not_oor= not_oor &  any( abs(fm_pts(:,3) ) <= z_depth , 2);
       end

       dims=1:2;

    elseif ec==4 %cm is 3D
       error('cant handle 3D coarse models (yet)'); 
    else
       error('coarse model must be 2 or 3D');
    end

    tsn(not_oor)= tsearchn(cm.nodes(:,dims), cm.elems, fm_pts(not_oor,dims));
    tsn= reshape( tsn, nf, df);


% interpolate over a triangle with n_interp points
% generate a set of points to fairly cover the triangle
% dim_coarse is dimensions + 1 of coarse model
function interp= triangle_interpolation(n_interp, dim_fine)
    interp= zeros(0,dim_fine);

    if dim_fine==3
       for i=0:n_interp
          for j=0:n_interp-i
             interp= [interp;i,j,n_interp-i-j];
          end
       end
    elseif dim_fine==4
       for i=0:n_interp
          for j=0:n_interp-i
             for k=0:n_interp-i-j
                interp= [interp;i,j,k,n_interp-i-j-k];
             end
          end
       end
    else
       error('cant handle dim_fine!=2');
    end

    interp= (interp + 1/dim_fine )/(n_interp+1);

% find fraction of elem contained in cm
% idx is the index into which the elem is contained
% if idx >= 1, the element is completely with in that coarse elem
% if idx == 0, the element is crosses several coarse elems
% if idx == Nan, the element is outside the coarse model
function c_elems = contained_elems_i( fm, cm, idx, z_depth)
   [nc,dc]= size(cm.elems);
   [nf,df]= size(fm.elems);

   fidx= find(idx==0);
   l_fidx= length(fidx);

   c_e_i= []; c_e_j=[]; c_e_v=[];

   % lower density interpolation in higher dimentions, since 
   % the added dimensions will give extra interpolation points.
   interp= triangle_interpolation( 7-df, df );

   l_interp = size(interp,1);
   dims = 1:dc-1; % run calc over dimensions 1 to dc-1

   el_nodes= fm.nodes(fm.elems(fidx,:)',:);
   % need to be of size df x N for reshape
   el_nodes= reshape(el_nodes, df, l_fidx*(df-1) );
   fm_pts = interp*el_nodes;
   fm_pts = reshape(fm_pts, l_fidx*l_interp, df-1);

   tsn= search_fm_pts_in_cm(cm, fm_pts, z_depth);

   tsn_idx= ones(l_interp,1)*fidx(:)';
   tsn_idx= tsn_idx(:)';
   % find and isolate Nans
   nan_idx= isnan(tsn);
   tsn(nan_idx) = [];
   tsn_idx(nan_idx) = [];
   % scale for effect of removed nans
   nan_weight= reshape(nan_idx, l_interp, l_fidx);
   nan_weight= l_interp - sum(nan_weight,1);
   
   in_idx= find((idx==0) | isnan(idx));
   ridx= 1:nf; ridx(in_idx)=[];
   idx(in_idx)=[];

   % first term is contribution from f_elems in one c_elem
   % next term is contribution from f_elems in many c_elems, weighted
   c_elems = sparse(ridx,idx,1,nf,nc) +  ...
             sparse(fidx,fidx,1./nan_weight,nf,nf) * ...
             sparse(tsn_idx,tsn,1,nf,nc);
   return
   c_elems1= sparse(ridx,idx,1,nf,nc);
   c_elems2= sparse(fidx,fidx,1./nan_weight,nf,nf);
   c_elems3= sparse(tsn_idx,tsn,1,nf,nc);
   

% non vectorized calculation
   for i = fidx'
      el_nodes= fm.nodes(fm.elems(i,:),:);
      fm_pts = interp*el_nodes;
      tsn= tsearchn(cm.nodes(:,dims), cm.elems, fm_pts(:,dims));
      tsn(isnan(tsn))=[];
      c_elems_i = sparse(1,tsn,1,1,nc)'/length(tsn);
      [c_i, c_j, c_v] = find(c_elems_i);
      c_e_i= [c_e_i;0*c_i+i];
      c_e_j= [c_e_j;c_j];
      c_e_v= [c_e_v;c_v];
%     c_elems(i,:)= sparse(1,tsn,1,1,nc)/length(tsn);
      if length(unique(tsn))==1;kk=kk+1; disp([kk,i]);end % how many unnecessary calcs?
   end

   ridx= 1:nf; ridx(fidx)=[];
   idx(fidx)=[];
   c_elems = sparse(ridx,idx,1,nf,nc) +  ...
             sparse(c_e_i,c_e_j,c_e_v,nf,nc);
      
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
