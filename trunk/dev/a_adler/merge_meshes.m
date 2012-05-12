function out = merge_meshes(M1,varargin)
% MERGE_MESHES - merges two meshes based on common nodes
% MERGE_MESHES(M1,M2,T) merges M2 in M1 using threshold T for detecting
%     corresponding nodes
% MERGE_MESHES(M1,M2,M3,..., T) merges M2, M3, ... into M1 (individually)

% (C) Bartlomiej Grychtol and Andy Adler, 2012. Licenced under GPL v2 or v3

th = varargin{end};
% Use a for loop, vectorised approach can run out of memory
l0 = length(M1.nodes);
if ~isfield(M1, 'mat_idx')
   M1.mat_idx = {1:length(M1.elems)};
end
for i = 1:length(varargin)-1
   l1 = length(M1.nodes);
   M2 = varargin{i};
   nodes_to_add = [];
   n_new_nodes = 0;
   match = 0 * (1:length(M2.nodes)); 
   for n = 1:length(M2.nodes)
      D = M1.nodes(1:l0,:) - repmat(M2.nodes(n,:),l0,1);
      D = sqrt(sum(D.^2,2));
      [val p] = min(D);
      if val < th
         match(n) = p;
      else
         match(n) = l1 + n_new_nodes+1;
         n_new_nodes = n_new_nodes + 1;
         nodes_to_add = [nodes_to_add; M2.nodes(n,:)];
      end
   end
   M1.nodes = [M1.nodes; nodes_to_add];
   M1.mat_idx = [M1.mat_idx {length(M1.elems)+(1:length(M2.elems))}];
   M1.elems = [M1.elems; match(M2.elems)];
   % this is not strictly correct, but visualizes nicely
   M1.boundary = [M1.boundary; match(M2.boundary)]; 
end

out =  M1;
% rmfield(M1,'boundary');
