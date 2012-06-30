function [mdl2,idx2] = mdl2d_from3d(mdl3,idx3);
   % set name
   mdl2 = eidors_obj('fwd_model',sprintf('%s 2D',mdl3.name));

   % set nodes
   [bdy,idx] = find_boundary(mdl3.elems);
   vtx = mdl3.nodes;
   z_vtx = reshape(vtx(bdy,3), size(bdy) );
   lay0  = find( all(z_vtx==0,2) );
   bdy0  = bdy( lay0, :);
   
   vtx0  = unique(bdy0(:));
   mdl2.nodes = vtx(vtx0,1:2);

   % set elems
   nmap  = zeros(size(vtx,1),1); nmap(vtx0) = 1:length(vtx0);
   bdy0  = reshape(nmap(bdy0), size(bdy0) ); % renumber to new scheme
   mdl2.elems = bdy0;

   % set boundary
   mdl2.boundary = find_boundary( mdl2.elems);

   % set gnd_node
   mdl2.gnd_node = nmap(mdl3.gnd_node);

   % set material indices
   % TODO: vectorize code
   idx2 = {};
   idx0  = idx( lay0, :);
   for i=1:size(idx3,2)
     idx2{i} = [];
     ii = 1;
     for j=1:size(idx3{i},1)
         idx_tmp = find( idx0==idx3{i}(j) );
         if not(isempty(idx_tmp))
           idx2{i}(ii,1) = idx_tmp(1,1);
           ii = ii + 1;
         end
     end
   end
   
   % set electrode
   if isfield(mdl3,'electrode')
     mdl2.electrode = mdl3.electrode;
     for i=1:length(mdl2.electrode);
        enodes = nmap( mdl2.electrode(i).nodes );
        enodes(enodes==0) = []; % Remove 3D layers
        mdl2.electrode(i).nodes = enodes(:)';
     end
   end

   % copy other fields
   if isfield(mdl3,'stimulation'); mdl2.stimulation= mdl3.stimulation; end
   %if isfield(mdl3,'solve');       mdl2.solve = mdl3.solve;            end
   mdl2.solve = 'aa_fwd_solve'; % FIXME? can't use default np_fwd_solve
   if isfield(mdl3,'jacobian');    mdl2.jacobian = mdl3.jacobian;      end
   %if isfield(mdl3,'system_mat');  mdl2.system_mat = mdl3.system_mat;  end
   mdl2.system_mat = 'aa_calc_system_mat'; % FIXME? can't use default np_calc_system_mat

   % update cache
   mdl2 = eidors_obj('fwd_model',mdl2);