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
% $Id$

mdl= eidors_obj('fwd_model', name);

for i= 1:prod(size(elec_nodes)) 
   vtx= [vtx; elec_nodes{i}];
end

% Sort and remove duplicates
vtx = unique(vtx, 'rows');

%vtx=  perturb_vtx(vtx) ;
%simp= delaunayn( vtx  );
simp= delaunayn( perturb_vtx(vtx)  );
% remove zero area elements
keep= ones(size(simp,1),1);
oo= ones(size(simp,2),1);
for i=1:size(simp,1);
   area= abs(det([oo,vtx(simp(i,:),:)]));
   if area<1e-6; keep(i)=0;end
end
simp= simp(find(keep),:);

% reorder elements so they're sorted by lowest node number
[jnk,idx] = sort(mean(simp,2));
simp= simp(idx,:);


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
mdl.solve=          @aa_fwd_solve;
mdl.jacobian=       @aa_calc_jacobian;
mdl.system_mat=     @aa_calc_system_mat;


function vtx_perturb= perturb_vtx( vtx );
if 0
   ctr= mean(vtx,1);
   v_ctr= vtx - ones(size(vtx,1),1) * ctr;
   r2 = (sum(v_ctr.^2,2)) .^ (1/4);
   vtx_perturb = vtx + .1* v_ctr ./ (r2*ones(1,size(vtx,2)));
else
   scl= std(vtx,1);
   vtx_perturb = vtx + 1e-5*randn(size(vtx)).*(ones(size(vtx,1),1)*scl);
end

   
function vtx_perturb= perturb_vtx_try1( vtx );
   max_vtx = max(vtx);
   min_vtx = min(vtx);
   vtx_perturb = vtx + 0.05*randn(size(vtx));
   for d= 1:size(vtx,2)
      idx= (vtx(:,d) == min_vtx(d)) | ...
           (vtx(:,d) == max_vtx(d));
      vtx_perturb(idx,d) = vtx(idx,d); 
   end
