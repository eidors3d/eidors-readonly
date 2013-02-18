function out = merge_meshes(M1,varargin)
%MERGE_MESHES - merges two meshes based on common nodes
% MERGE_MESHES(M1,M2,T) merges M2 in M1 using threshold T for detecting
%     corresponding nodes. The meshes must not overlap.
% MERGE_MESHES(M1,M2,M3,..., T) merges M2, M3, ... into M1 (individually)
%
% Note that the boundaries of the separate meshes will only be
% concatenated, as this visualises nicely. To calculate the correct
% boundary use FIND_BOUNDARY.
%
% See also FIND_BOUNDARY

% (C) Bartlomiej Grychtol and Andy Adler, 2012-2013. Licence: GPL v2 or v3
% $Id$

if ischar(M1) && strcmp(M1,'UNIT_TEST'), run_unit_test; return; end

if nargin < 3  || isstruct(varargin{end})
   th = mean(std(M1.nodes))/length(M1.nodes);
   shapes = varargin;
else
   th = varargin{end};
   shapes = varargin(1:end-1);
end
% Use a for loop, vectorised approach can run out of memory
l0 = length(M1.nodes);
if ~isfield(M1, 'mat_idx')
   M1.mat_idx = {1:length(M1.elems)};
end
for i = 1:length(shapes)
   l1 = length(M1.nodes);
   M2 = shapes{i};
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

function run_unit_test
subplot(221)
cyl = ng_mk_cyl_models(3,[0],[]);
show_fem(cyl)

subplot(222)
top_nodes = cyl.nodes(:,3)>=1.5;
top_elems = sum(top_nodes(cyl.elems),2)==4;
top.elems = cyl.elems(top_elems,:);
nds = unique(top.elems);
map = zeros(1,length(cyl.nodes));
map(nds) = 1:length(nds);
top.elems = map(top.elems);
top.nodes = cyl.nodes(nds,:);
top.type = 'fwd_model';
top.boundary = find_boundary(top);
show_fem(top)
zlim([0 3]);

subplot(223)
bot_elems = ~top_elems;
bot.elems = cyl.elems(bot_elems,:);
nds = unique(bot.elems);
map = zeros(1,length(cyl.nodes));
map(nds) = 1:length(nds);
bot.elems = map(bot.elems);
bot.nodes = cyl.nodes(nds,:);
bot.type = 'fwd_model';
bot.boundary = find_boundary(bot);
show_fem(bot)
zlim([0 3]);


subplot(224)
M = merge_meshes(bot, top);
show_fem(M);

unit_test_cmp('Number of nodes',length(cyl.nodes), length(M.nodes),0);
unit_test_cmp('Number of elems',length(cyl.elems), length(M.elems),0);
