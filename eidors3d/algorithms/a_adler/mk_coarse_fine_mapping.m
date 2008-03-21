function [mapping,c_elems] = mk_coarse_fine_mapping( f_mdl, c_mdl );
% MK_COARSE_FINE_MAPPING: create a mapping matrix from coarse to fine FEM
% mapping= mk_coarse_fine_mapping( f_mdl, c_mdl );
%
% Mapping approximates elem_data_fine from elem_data_coase
%   elem_data_fine = Mapping*elem_data_coase
%
% c_mdl is coarse fwd_model
% f_mdl is fine fwd_model
% f_mdl.mk_coarse_fine_mapping.n_interp - default 50
%                 - number of points to interpolate in each dimension

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_coarse_fine_mapping.m,v 1.13 2008-03-21 01:24:03 aadler Exp $

% Mapping depends f_mdl and c_mdl, but only on nodes and elems
cc_mdl.elems = c_mdl.elems;
cto3d = speye(size(c_mdl.nodes,2), 3);
cc_mdl.nodes = c_mdl.nodes * cto3d;
cc_mdl.type  = 'fwd_model';
n_cc= size(c_mdl.elems, 1);

ff_mdl.elems = f_mdl.elems;
fto3d = speye(size(f_mdl.nodes,2), 3);
ff_mdl.nodes = f_mdl.nodes * fto3d;
ff_mdl.type  = 'fwd_model';
n_ff= size(f_mdl.elems, 1);

mapping = eidors_obj('get-cache', {ff_mdl,cc_mdl}, 'coarse_fine_mapping');
if ~isempty(mapping)
    eidors_msg('mk_coarse_fine_mapping: using cached value', 3);
else

    try
       n_interp = f_mdl.mk_coarse_fine_mapping.n_interp;
    catch
       n_interp = 50;
    end
   

    eidors_obj('set-cache', {ff_mdl,cc_mdl}, 'coarse_fine_mapping', mapping);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end

f_elems = all_contained_elems( ff_mdl, cc_mdl);
c_elems = contained_elems_i( ff_mdl, cc_mdl, f_elems);

function mapping= mk_mapping( ff_mdl, cc_mdl, n_interp);
    xyzmin = min([ff_mdl.nodes;cc_mdl.nodes]);
    xyzmax = max([ff_mdl.nodes;cc_mdl.nodes]);
    xyz = interpxyz( xyzmin, xyzmax, n_interp);

    f_tria= mk_tri_pts( ff_mdl, xyz, fto3d);
    c_tria= mk_tri_pts( cc_mdl, xyz, cto3d);

    ff= find(~isnan(f_tria) & ~isnan(c_tria));
    c2f = sparse(f_tria(ff),c_tria(ff),1,n_ff, n_cc);
    csum = sparse(f_tria(ff),1,1,n_ff, 1);
    mapping = c2f./(csum*ones(1,n_cc));

% find all elems of ff_mdl completely contained in cc_mdl
function c_elems = all_contained_elems( fm, cm)
    nd= size(fm.nodes,2); % n dims
    [nf,ef]= size(fm.elems);
    nc= size(cm.elems,1);
    fm_pts = zeros(nf*ef,3);
    % shrink pts slightly so they're not on the boundary
    s_fac= 1-1e-6;
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
    tsn= sort( tsn, 2); % send Nans to right
    isn= find(isnan(tsn));
    tsn(isn)= tsn( 1+rem(isn-1,nf ) );
    c_elems= all(diff(tsn,1,2)==0,2);
    c_elems= c_elems.* tsn(:,1);

% interpolate over a triangle with n_interp points
% generate a set of points to fairly cover the triangle
function interp= triangle_interpolation(n_interp)
    interp= zeros(0,3);
    for i=0:n_interp
       for j=0:n_interp-i
          interp= [interp;i,j,n_interp-i-j];
       end
    end

    interp= (interp+1/3)/(n_interp+1);

% find fraction of elem contained in cm
function c_elems = contained_elems_i( fm, cm, idx)
   nc= size(cm.elems,1);
   nf= size(fm.elems,1);

   fidx= find(idx==0);
   ridx= 1:nf; ridx(fidx)=[];
   idx(fidx)=[];
   c_elems = sparse(ridx,idx,1,nf,nc);

   interp= triangle_interpolation( 4 );
   l_interp= size(interp,1);
   for i = fidx'
      el_nodes= fm.nodes(fm.elems(i,:),:);
      fm_pts = interp*el_nodes;
      tsn= tsearchn(cm.nodes(:,1:2), cm.elems, fm_pts(:,1:2));
      c_elems(i,:)= sparse(1,tsn,1,1,nc)/l_interp;
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


% calculate mapping of points xyz into triangle
% fmdl is fwd_mode, xyz is points [x(:), y(:), z(:)]
% to 3d is 
function tri_pts = mk_tri_pts( fmdl, xyz, to3d )
    nodes = fmdl.nodes*to3d'; 
    elems = fmdl.elems;
    tri_pts= tsearchn(nodes, elems, xyz*to3d');

    ff= find(~isnan(tri_pts));
    n_ff= size(elems, 1);
    csum = sparse(tri_pts(ff),1,1,n_ff, 1);
    if any(csum<10)
       warning(sprintf( ...
         'mesh density is too low: min=%d', full(min(csum)) ))
    end
