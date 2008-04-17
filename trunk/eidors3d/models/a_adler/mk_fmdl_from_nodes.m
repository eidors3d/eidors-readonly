function mdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, name)
% MK_FMDL_FROM_NODES: create fmdl from nodes
% fwd_mdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, name)
%
% Create a fwd_model by delaynay triangularization of vtx nodes
%  vtx:         matrix of nodes in model (can (or not) include elec_nodes)
%  elec_nodes:  cell of matrix N x [x,y,{z}] for each electrode
%  z_contact:   vector or scalar electrode contact impedance
%  name:        name for eidors object
%
%  fwd_mdl:     eidors format fwd_model
%
% Limitations: assumes a convex model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_fmdl_from_nodes.m,v 1.1 2008-04-17 15:54:56 aadler Exp $

mdl= eidors_obj('fwd_model', name);

for i= 1:prod(size(elec_nodes)) 
   vtx= [vtx; elec_nodes{i}];
end

% Sort and remove duplicates
vtx = unique(vtx, 'rows');
simp= delaunayn( vtx );


mdl.nodes    = vtx;
mdl.elems    = simp;
mdl.boundary = find_boundary( simp );
mdl.gnd_node = 1;
mdl.np_fwd_solve.perm_sym = '{n}';

% Electrodes and z_contact

n_elec= prod(size(elec_nodes));
z_contact= z_contact.*ones(n_elec,1);
for i= 1:n_elec
   [jnk, idxa]            = intersect( vtx, elec_nodes{i},'rows');
   electrodes(i).nodes    = idxa(:)';
   electrodes(i).z_contact= z_contact(i);
end


mdl.electrode =     electrodes;
mdl.solve=          'np_fwd_solve';
mdl.jacobian=       'np_calc_jacobian';
mdl.system_mat=     'np_calc_system_mat';
