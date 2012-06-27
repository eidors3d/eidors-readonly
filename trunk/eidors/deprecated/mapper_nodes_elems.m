function EtoN= mapper_elems_nodes(fwd_model)
% MAPPER_ELEMS_NODES: calculates mapping function taking elems to nodes
% EtoN= mapper_elems_nodes(fwd_model)
% EtoN is a matrix [Nnodes x Nelems]
%    where [EtoN]_ij = 1/3 if elem_j shares node_i

warning('EIDORS:deprecated','MAPPER_ELEMS_NODES is deprecated as of 27-Jun-2012.');

elems= fwd_model.elems;
n_elems= size(elems,1);
n_dimp1= size(elems,2);
n_nodes= size(fwd_model.nodes,1);
EtoN = sparse(elems, (1:n_elems)'*ones(1,n_dimp1), 1/3, n_nodes, n_elems);
