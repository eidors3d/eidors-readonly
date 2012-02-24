function [bound,elem,nodes]=mc_fem_modify(fwd_model)
%Function  that takes a fwd_model structure (with at least fwd_model.boundary,
%fwd_model.elems, fwd_model.nodes), reorders these to form anti-clockwise sets, and
%then adds extra nodes for p-refinement.
%
%Function exits with fwd_model.elem(i).nodes; fwd_model.bound(i).nodes, which have the
%extra p-refined nodes correctly positioned, and still has the fwd_model.boundary
%and fwd_model.elems structures

%Change to fwd_model.elem(i)
for i=1:size(fwd_model.elems,1)
    fwd_model.elem(i).nodes=fwd_model.elems(i,:);
end

%Change to fwd_model.boundary(i)
for i=1:size(fwd_model.boundary,1)
    fwd_model.bound(i).nodes=fwd_model.boundary(i,:);
end

%Prefine - boundary first, then remaining elements
fwd_model=mc_p_refine(fwd_model);

%FIXME Clear this up, by changing the output of two functions above, 
%so that the output is just bound and elem...
bound=fwd_model.bound; elem=fwd_model.elem; nodes=fwd_model.nodes;