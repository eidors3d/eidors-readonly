function mapping= mk_coarse_fine_mapping( f_mdl, c_mdl );
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

% (C) 2007 Andy Adler. Licenced under the GPL Version 2
% $Id: mk_coarse_fine_mapping.m,v 1.8 2007-08-29 09:26:17 aadler Exp $

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

mapping = eidors_obj('get-cache', ff_mdl, 'coarse_fine_mapping', cc_mdl);
if ~isempty(mapping)
    eidors_msg('mk_coarse_fine_mapping: using cached value', 3);
else

    try
       n_interp = f_mdl.mk_coarse_fine_mapping.n_interp;
    catch
       n_interp = 50;
    end
   
    xyzmin = min([ff_mdl.nodes;cc_mdl.nodes]);
    xyzmax = max([ff_mdl.nodes;cc_mdl.nodes]);
    xyz = interpxyz( xyzmin, xyzmax, n_interp);

    f_tria= mk_tri_pts( ff_mdl, xyz, fto3d);
    c_tria= mk_tri_pts( cc_mdl, xyz, cto3d);

    ff= find(~isnan(f_tria));
    c2f = sparse(f_tria(ff),c_tria(ff),1,n_ff, n_cc);
    csum = sparse(f_tria(ff),1,1,n_ff, 1);
    mapping = c2f./(csum*ones(1,n_cc));

    eidors_obj('set-cache', ff_mdl, 'coarse_fine_mapping', mapping, cc_mdl);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
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
