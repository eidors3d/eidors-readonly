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

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id$

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
    xyzdelta= xyzmax - xyzmin;
    xyz_interp = 1 + floor(n_interp * xyzdelta / max(xyzdelta) );
    xspace = linspace(xyzmin(1), xyzmax(1), xyz_interp(1) );
    yspace = linspace(xyzmin(2), xyzmax(2), xyz_interp(2) );
    zspace = linspace(xyzmin(3), xyzmax(3), xyz_interp(3) );
    [xx3,yy3,zz3] = ndgrid( xspace, yspace, zspace );
    xyz= [xx3(:), yy3(:), zz3(:)];

    f_tria= tsearchn(ff_mdl.nodes*fto3d', ff_mdl.elems, xyz*fto3d');
    c_tria= tsearchn(cc_mdl.nodes*cto3d', cc_mdl.elems, xyz*cto3d');

    ff= find(~isnan(f_tria));
    c2f = sparse(f_tria(ff),c_tria(ff),1,n_ff, n_cc);
    csum = sparse(f_tria(ff),1,1,n_ff, 1);
    mapping = c2f./(csum*ones(1,n_cc));

    eidors_obj('set-cache', ff_mdl, 'coarse_fine_mapping', mapping, cc_mdl);
    eidors_msg('mk_coarse_fine_mapping: setting cached value', 3);
end
