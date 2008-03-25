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
% indicates the +/- z_depth which is reflected onto the
% c_mdl (ie c_mdl is at z=0 +/- z_depth)

% (C) 2007-2008 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_coarse_fine_mapping.m,v 1.17 2008-03-25 14:37:46 aadler Exp $

% Mapping depends only on nodes and elems
try; c_mdl= rmfield(c_mdl,'electrode');   end
try; c_mdl= rmfield(c_mdl,'stimulation'); end
try; f_mdl= rmfield(f_mdl,'electrode');   end
try; f_mdl= rmfield(f_mdl,'stimulation'); end

% TODO: mapping step-> map f_mdl onto c_mdl

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

    % look for f_mdl z out of range
    [nf,ne] = size(f_mdl.elems);
    f_elem_z= reshape(f_mdl.nodes( f_mdl.elems(:), 3), nf, ne);
    oor= all( abs(f_elem_z) > z_depth , 2);

    mapping = spdiags(~oor, 0, nf,nf) * mapping;

    eidors_obj('set-cache', {f_mdl,c_mdl}, 'coarse_fine_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end


% find all elems of ff_mdl completely contained in cc_mdl
function c_elems = all_contained_elems( fm, cm, z_depth)
    nd= size(fm.nodes,2); % n dims
    [nf,ef]= size(fm.elems);
    nc= size(cm.elems,1);
    fm_pts = zeros(nf*ef,3);
    % shrink pts slightly so they're not on the boundary
    s_fac= .9999;
    for dim= 1:nd
       % fm_pts is local_nodes x elems x xyz
       fm_pt= reshape(fm.nodes(fm.elems,dim),nf,ef);
       fm_ctr= mean(fm_pt,2)*ones(1,ef);
       fm_pt = s_fac*(fm_pt-fm_ctr) + fm_ctr;
       fm_pts(:,dim) = fm_pt(:);
    end

    tsn= tsearchn(cm.nodes(:,1:2), cm.elems, fm_pts(:,1:2));
    tsn= reshape( tsn, nf, ef);
    % if node in fm is outside cm, then set it back inside
    % this way only elems in fm that are completely outside cm
    % will cause a problem
    tsn= sort( tsn, 2); % send Nans to right
    isn= find(isnan(tsn));
    tsn(isn)= tsn( 1+rem(isn-1,nf ) );
    c_elems= all(diff(tsn,1,2)==0,2);
    c_elems= c_elems.* tsn(:,1);

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
function c_elems = contained_elems_i( fm, cm, idx)
   [nc,dc]= size(cm.elems);
   [nf,df]= size(fm.elems);

   fidx= find(idx==0);
   ridx= 1:nf; ridx(fidx)=[];
   idx(fidx)=[];
   c_elems = sparse(ridx,idx,1,nf,nc);

   if df==3
      interp= triangle_interpolation( 4, df );
   else
      interp= triangle_interpolation( 3, df );
   end
   dims = 1:dc-1; % run calc over dimensions 1 to dc-1
   for i = fidx'
      el_nodes= fm.nodes(fm.elems(i,:),:);
      fm_pts = interp*el_nodes;
      tsn= tsearchn(cm.nodes(:,dims), cm.elems, fm_pts(:,dims));
      tsn(isnan(tsn))=[];
      c_elems(i,:)= sparse(1,tsn,1,1,nc)/length(tsn);
%     if length(unique(tsn))==1; disp(i);end % how many unnecessary calcs?
   end
      
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
