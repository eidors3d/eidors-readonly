function fmdl = remove_unused_nodes( fmdl );
% REMOVE_UNUSED_NODES: identify and remove unused nodes in model
% Usage: fmdl = remove_unused_nodes( fmdl );

% (C) 2019 Andy Adler. License: GPL v2 or v3. $Id$

   if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

   if num_elems(fmdl)==0; return; end; % don't operate on pathalogical models

   usednodes = unique(fmdl.elems(:));
   if max(usednodes) > num_nodes(fmdl)
      error('remove_unused_nodes: more nodes are used than exist');
   end
   nidx = zeros(num_nodes(fmdl),1);
   nidx(usednodes) = 1:length(usednodes);
   fmdl.nodes(nidx==0,:) = [];
   fmdl.elems = remap(nidx, fmdl.elems);

   for i=1:length(fmdl.electrode)
      fmdl.electrode(i).nodes =  remap(nidx, fmdl.electrode(i).nodes);
      removed = fmdl.electrode(i).nodes == 0;
      fmdl.electrode(i).nodes( removed ) = [];
      if isfield(fmdl.electrode(i),'faces')
          fmdl.electrode(i).faces =  remap(nidx, fmdl.electrode(i).faces);
          removed = any(fmdl.electrode(i).face == 0,2);
          fmdl.electrode(i).faces(removed,:) = [];
          if isempty(fmdl.electrode(i).faces);
             eidors_msg('Zeros in faces #%d',i,1);
             keyboard
          end
      else
          if isempty(fmdl.electrode(i).nodes);
             eidors_msg('Zeros in nodes #%d',i,1);
             keyboard
          end
      end
   end
%  fmdl.boundary = find_boundary(fmdl);
   fmdl.boundary = remap(nidx, fmdl.boundary);
   fmdl.boundary(any(fmdl.boundary==0,2),:) = [];
   fmdl.gnd_node = nidx(fmdl.gnd_node);
   if fmdl.gnd_node == 0 %% New gnd node if missing
      fmdl = assign_new_gnd_node( fmdl );
   end

%  fmdl.elems = reshape(nidx(fmdl.elems),[],size(fmdl.elems,2));
function mat = remap(nidx,mat)
   mat = reshape(nidx(mat),size(mat));

function fmdl = assign_new_gnd_node( fmdl );
   eidors_msg('FEM_ELECTRODE: Lost ground node => replacing',1);
   d2 = sum((fmdl.nodes - repmat(mean(fmdl.nodes, 1), size(fmdl.nodes, 1), 1)).^2,2);
   [~,fmdl.gnd_node] = min(d2);

function do_unit_test
   fmdl = getfield(mk_common_model('a2c2',4),'fwd_model');
   fmdl.elems(1:4,:) = [];
   fmdl = remove_unused_nodes(fmdl);

   unit_test_cmp('nodes',size(fmdl.nodes),[40,2]);
   unit_test_cmp('elems',size(fmdl.elems),[60,3]);
   unit_test_cmp('ground',fmdl.gnd_node,3);
