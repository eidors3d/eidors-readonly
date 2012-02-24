function mdl=mc_fem_modify(mdl,eletype)
%Function  that takes a fwd_model structure (with at least mdl.boundary,
%mdl.elems, mdl.nodes), reorders these to form anti-clockwise sets, and
%then adds extra nodes for p-refinement.
%
%Function exits with mdl.elem(i).nodes; mdl.bound(i).nodes, which have the
%extra p-refined nodes correctly positioned, and still has the mdl.boundary
%and mdl.elems structures

%Give the model a finite element type
mdl.mc_type = eletype;

%Change to mdl.elem(i)
for i=1:size(mdl.elems,1)
    mdl.elem(i).nodes=mdl.elems(i,:);
end

%Change to mdl.boundary(i)
for i=1:size(mdl.boundary,1)
    mdl.bound(i).nodes=mdl.boundary(i,:);
end

%Reorder nodes in element 
mdl=mc_elem_bound_reorder(mdl);

%Prefine - boundary first, then remaining elements
mdl=mc_p_refine(mdl);