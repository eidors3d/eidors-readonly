function c2f = mk_tri2tet_c2f(fmdl,rmdl, opt)
%MK_TRI2TET_C2F - coarse2fine mapping between tri-based and tet-based models
% C2F = MK_TRI2TET_C2F(FMDL,RMDL) returns in C2F the fraction of volume of
% each element of the fine (tet-based) model contained in each element of
% the coarse (tri-based) model. Each triangle is extruded along the z-axis
% to form a prism.
% Uses CONVHULLN to calculate the volume defined by a set of intersection
% points between individual tet and vox elements.
%
% C2F = MK_TRI2TET_C2F(FMDL,RMDL,OPT) allows specifying options.
% 
% Inputs:
%   FMDL - an EIDORS (tet-based) 3D forward model
%   RMDL - an EIDORS (tri-based) 2D forward model
%   OPT  - an option structure with the following fields and defaults:
%      .z_depth
%      .do_not_scale  - set to true to prevent scaling the models to unit
%                       cube before any calculations, including thresholds.
%                       Default: false
%      .tol_node2tri  - minimum value of a barycentric coordinate to 
%                       decide a tet node is lying inside a triangle and 
%                       not on its edge. Default: eps
%      .tol_node2tet  - tolerance for determinant <= 0 in testing for
%                       points inside tets. Default: eps
%      .tol_edge2edge - maximum distance between "intersecting" edges
%                       Default: 6*sqrt(3)*eps(a), where a is
%                       min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))
%      .tol_edge2tri  - minimum value of a barycentric coordinate to 
%                       decide an edge-plane intersection point is lying 
%                       inside a triangle. Default: eps
%                       
% NOTE that for grid-based models, such as returned by MK_GRID_MODEL or
% MK_VOXEL_VOLUME, MK_GRID_C2F is much faster.
%
% Set eidors_msg 'log level' < 2 to supress output to command line.
% 
% See also MK_GRID_C2F, MK_TET_C2F, MK_TRI_C2F, MK_COARSE_FINE_MAPPING,
%          FIND_EDGE2EDGE_INTERSECTIONS, CONVHULLN, MK_APPROX_C2F, 
%          POINT_IN_TRIANGLE, EIDORS_MSG

% (C) 2015 Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$


if nargin == 0 || (ischar(fmdl) && strcmp(fmdl,'UNIT_TEST')) 
   do_unit_test; 
   return; 
end
if nargin < 3
   opt = struct();
end

f_elems = size(fmdl.elems,1);
r_elems = size(rmdl.elems,1);
c2f = sparse(f_elems,r_elems);

if size(rmdl.nodes,2) == 2
   rmdl.nodes(:,3) = max(fmdl.nodes(:,3))/2 + min(fmdl.nodes(:,3))/2;
end

[fmdl,rmdl,fmdl_idx,rmdl_idx] = crop_models(fmdl,rmdl, opt);

if ~(any(fmdl_idx) && any(rmdl_idx))
   eidors_msg('@@: models do not overlap, returning all-zeros');
   return
end

opt = parse_opts(fmdl,rmdl, opt);

[fmdl,rmdl, opt] = center_scale_models(fmdl,rmdl, opt);

copt.fstr = 'mk_tri2tet_c2f';
c2f(fmdl_idx,rmdl_idx) = eidors_cache(@do_mk_tri2tet_c2f,{fmdl,rmdl,opt},copt);
% c2f = eidors_cache(@do_mk_tri2tet_c2f,{fmdl,rmdl,opt},copt);

end

function c2f = do_mk_tri2tet_c2f(fmdl,rmdl,opt)

   DEBUG =  eidors_debug('query','mk_tri2tet_c2f');

   r_elems = size(rmdl.elems,1);
   f_elems = size(fmdl.elems,1);
 
   c2f = sparse(f_elems, r_elems);
   
   progress_msg('Prepare tet model...');
   fmdl = prepare_tet_mdl(fmdl);
   progress_msg(Inf);
   
   progress_msg('Prepare tri model...');
   rmdl = prepare_tri_mdl(rmdl);
   progress_msg(Inf);

   pmopt.final_msg = 'none';
   if ~isinf(opt.z_depth)
      % tri nodes on top and bot planes inside tets
      progress_msg('Find top cap nodes in tets...',-1,pmopt);
      [top_node2tet, top_nodes, top_nodes_above, top_nodes_below] = ...
         get_cap_nodes_in_tets(fmdl,rmdl,opt.top,opt.tol_node2tet);
      progress_msg(sprintf('Found %d', nnz(top_node2tet)), Inf);
      
      progress_msg('Find bottom cap nodes in tets...',-1,pmopt);
      [bot_node2tet, bot_nodes, bot_nodes_above, bot_nodes_below] = ...
         get_cap_nodes_in_tets(fmdl,rmdl,opt.bot,opt.tol_node2tet);
      progress_msg(sprintf('Found %d', nnz(bot_node2tet)), Inf);
      
      
      progress_msg('Find tet_face2tri_edge intersections (top) ...');
      [intpts4, top_tet_face2tri_edge,top_tet_face2intpt4,top_tri_edge2intpt4] = ...
         get_cap_tet_face2tri_edge_intersections( fmdl,rmdl,opt.top, ...
         top_nodes_above,top_nodes_below, opt.tol_edge2tri);
      progress_msg(sprintf('Found %d', size(intpts4,1)),Inf);
      
      progress_msg('Find tet_face2tri_edge intersections (bot) ...');
      [intpts5, bot_tet_face2tri_edge,bot_tet_face2intpt5,bot_tri_edge2intpt5] = ...
         get_cap_tet_face2tri_edge_intersections( fmdl,rmdl,opt.bot, ...
         bot_nodes_above,bot_nodes_below, opt.tol_edge2tri);
      progress_msg(sprintf('Found %d', size(intpts5,1)),Inf);
      
      % find tet_edge2tri_face intersections
      progress_msg('Find tet_edge2tri_face intersections (top) ...');
      [intpts6, top_tet_edge2tri_face, top_tet_edge2intpt6, tri2intpt6] = ...
         get_cap_tet_edge2tri_face_intersections(fmdl,rmdl,opt.top, ...
         top_nodes_above, top_nodes_below, opt.tol_edge2tri);
      progress_msg(sprintf('Found %d', size(intpts6,1)),Inf);
      
      
      progress_msg('Find tet_edge2tri_face intersections (bot) ...');
      [intpts7, bot_tet_edge2tri_face, bot_tet_edge2intpt7, tri2intpt7] = ...
         get_cap_tet_edge2tri_face_intersections(fmdl,rmdl,opt.bot, ...
         bot_nodes_above, bot_nodes_below, opt.tol_edge2tri);
      progress_msg(sprintf('Found %d', size(intpts6,1)),Inf);
      
      
      
      progress_msg('Find tet_edge2tri_edge intersections (top) ...',-1,pmopt);
      edgeidx = any(top_nodes_above(fmdl.edges),2) & any(top_nodes_below(fmdl.edges),2);
      nodes = rmdl.nodes;
      nodes(:,3) = opt.top;
      top_tet_edge2tri_edge = sparse(size(fmdl.edges,1),size(rmdl.edges,1));
      [intpts8, top_tet_edge2tri_edge(edgeidx,:), top_tet_edge2intpt8(edgeidx,:), top_tri_edge2intpt8] = ...
         find_edge2edge_intersections(fmdl.edges(edgeidx,:),fmdl.nodes, ...
         rmdl.edges,nodes, opt.tol_edge2edge);
      if size(top_tet_edge2intpt8,1)~=size(fmdl.edges,1)
         if ~isempty(intpts8)
            top_tet_edge2intpt8(size(fmdl.edges,1),1) = 0;
         else
            top_tet_edge2intpt8 = sparse(size(fmdl.edges,1),0);
         end
      end
      progress_msg(sprintf('Found %d', size(intpts8,1)),Inf);
      
      
      progress_msg('Find tet_edge2tri_edge intersections (bot) ...',-1,pmopt);
      edgeidx = any(bot_nodes_above(fmdl.edges),2) & any(bot_nodes_below(fmdl.edges),2);
      nodes = rmdl.nodes;
      nodes(:,3) = opt.bot;
      bot_tet_edge2tri_edge = sparse(size(fmdl.edges,1),size(rmdl.edges,1));
      [intpts9, bot_tet_edge2tri_edge(edgeidx,:), ...
         bot_tet_edge2intpt9(edgeidx,:), bot_tri_edge2intpt9] = ...
         find_edge2edge_intersections(fmdl.edges(edgeidx,:),fmdl.nodes, ...
         rmdl.edges,nodes, opt.tol_edge2edge);
      if size(bot_tet_edge2intpt9,1)~=size(fmdl.edges,1)
         if ~isempty(intpts9)
            bot_tet_edge2intpt9(size(fmdl.edges,1),1) = 0;
         else
            bot_tet_edge2intpt9 = sparse(size(fmdl.edges,1),0);
         end
      end
      progress_msg(sprintf('Found %d', size(intpts9,1)),Inf);
      
      % tri contained in tet
      progress_msg('Find tris contained in tet...')
      tri_in_tet = rmdl.node2elem' * (bot_node2tet + top_node2tet) == 6;
      progress_msg(sprintf('Found %d',nnz(tri_in_tet)), Inf);
   end

   % tet nodes inside triangles
   % positive tolerance to catch also points on the edges.
   progress_msg('Find tet_nodes in tri_elems...');
   in = (fmdl.nodes(:,3) >= opt.bot) & (fmdl.nodes(:,3) <= opt.top);
   tet_node2tri = spalloc(size(fmdl.nodes,1),size(rmdl.elems,1),nnz(in));
   tet_node2tri(in,:) = point_in_triangle(fmdl.nodes(in,1:2), ...
                                    rmdl.elems, ...
                                    rmdl.nodes(:,1:2), ...
                                    opt.tol_node2tri);      
   progress_msg(sprintf('Found %d', nnz(tet_node2tri)), Inf);
  
   n_nodes = size(rmdl.nodes,1);
   vert_edges(:,1) = 1:n_nodes;
   vert_edges(:,2) = vert_edges(:,1) + n_nodes;
   
   progress_msg('Find tet_edge2vert_edge intersections...',-1,pmopt);
   if isinf(opt.z_depth)
      bot_nodes = rmdl.nodes; bot_nodes(:,3) = opt.bot - 1;
      top_nodes = rmdl.nodes; top_nodes(:,3) = opt.top + 1;
   end
   [intpts1,tet_edge2vert_edge,tet_edge2intpt1,vert_edge2intpt1] = ...
      find_edge2edge_intersections( fmdl.edges, fmdl.nodes,...
                                    vert_edges, [bot_nodes; top_nodes], ...
                                    opt.tol_edge2edge);
   progress_msg(sprintf('Found %d',size(intpts1,1)),Inf);
   
   % find tet_edges that cross the ractangular faces between the top and
   % bottom caps
   progress_msg('Find tet_edge2vert_face intersections...');
   if isinf(opt.z_depth)
      top = opt.top + 1;
      bot = opt.bot - 1;
   else
      top = opt.top;
      bot = opt.bot;
   end
   [intpts2,tet_edge2tri_edge,tet_edge2intpt2,tri_edge2intpt2] = ...
      find_edge2edge_intersections_2d( fmdl.edges, fmdl.nodes(:,1:3), ...
                                       rmdl.edges, rmdl.nodes(:,1:2), ...
                                       top, bot);
   progress_msg(sprintf('Found %d',size(intpts2,1)),Inf);

   
   progress_msg('Find tet_face2vert_edge intersections...');
   [intpts3, tet_face2vert_edge, tet_face2intpt3, vert_edge2intpt3] = ...
      get_tet_intersection_pts(fmdl,rmdl,top,bot,opt.tol_edge2tri);
   progress_msg(sprintf('Found %d',size(intpts3,1)),Inf);
   
      
   % tet contained in tri
   progress_msg('Find tets contained in tri...');
   tet_in_tri = (double(fmdl.node2elem') * tet_node2tri) == 4;
   progress_msg(sprintf('Found %d',nnz(tet_in_tri)), Inf);
   
   % tets and tris that intersect
   progress_msg('Find total tri v. tet intersections...');
   tri2intTet = double(rmdl.edge2elem') * tet_edge2tri_edge' * fmdl.edge2elem ...
      |double(rmdl.node2elem') * (tet_face2vert_edge>0)' * fmdl.elem2face';
   if ~isinf(opt.z_depth)
      tri2intTet = tri2intTet ...
      | double(rmdl.edge2elem') * (top_tet_face2tri_edge + bot_tet_face2tri_edge)' * fmdl.elem2face' ...
      | (top_tet_edge2tri_face + bot_tet_edge2tri_face)' * fmdl.edge2elem;
   end
   % exclude complete inclusion
   tri2intTet = tri2intTet & ~tet_in_tri';
   if ~isinf(opt.z_depth)
      tri2intTet = tri2intTet & ~tri_in_tet;
   end
   progress_msg(sprintf('Found %d',nnz(tri2intTet)), Inf);
   
   
   progress_msg('Calculate intersection volumes...');
   % sparse logical multiplication doesn't exist
   tri2intpt1 = logical(rmdl.node2elem'*vert_edge2intpt1)';
   tet2intpt1 = logical(fmdl.edge2elem'* tet_edge2intpt1)';
   
   tet2intpt2 = logical(fmdl.edge2elem'*tet_edge2intpt2)';
   tri2intpt2 = logical(rmdl.edge2elem'*tri_edge2intpt2)';
   
   tet2intpt3 = logical(double(fmdl.elem2face)*tet_face2intpt3)';
   tri2intpt3 = logical(rmdl.node2elem'*vert_edge2intpt3)';
   if ~isinf(opt.z_depth)
      tet2intpt4  = logical(double(fmdl.elem2face)*top_tet_face2intpt4)';
      tri2intpt4  = logical(rmdl.edge2elem'*top_tri_edge2intpt4)';
      
      tet2intpt5  = logical(double(fmdl.elem2face)*bot_tet_face2intpt5)';
      tri2intpt5  = logical(rmdl.edge2elem'*bot_tri_edge2intpt5)';
      
      tet2intpt6  = logical(fmdl.edge2elem'*top_tet_edge2intpt6)';
      tri2intpt6  = tri2intpt6';
      
      tet2intpt7  = logical(fmdl.edge2elem'*bot_tet_edge2intpt7)';
      tri2intpt7  = tri2intpt7';
   
      tet2intpt8  = logical(fmdl.edge2elem'*top_tet_edge2intpt8)';
      tri2intpt8  = logical(rmdl.edge2elem'*top_tri_edge2intpt8)';
      
      tet2intpt9  = logical(fmdl.edge2elem'*bot_tet_edge2intpt9)';
      tri2intpt9  = logical(rmdl.edge2elem'*bot_tri_edge2intpt9)';
   end
   tri_todo = find(sum(tri2intTet,2)>0);
   C = []; F = []; V = [];
   
   id = 0; lvox = length(tri_todo);
   mint = ceil(lvox/100);
   
   for v = tri_todo'
      id = id+1;
      if mod(id,mint)==0, progress_msg(id/lvox); end
      tet_todo = find(tri2intTet(v,:));
      common_intpts1 = bsxfun(@and,tri2intpt1(:,v), tet2intpt1(:,tet_todo));
      common_intpts2 = bsxfun(@and,tri2intpt2(:,v), tet2intpt2(:,tet_todo));
      common_intpts3 = bsxfun(@and,tri2intpt3(:,v), tet2intpt3(:,tet_todo));
      if ~isinf(opt.z_depth)
         common_intpts4 = bsxfun(@and,tri2intpt4(:,v), tet2intpt4(:,tet_todo));
         common_intpts5 = bsxfun(@and,tri2intpt5(:,v), tet2intpt5(:,tet_todo));
         common_intpts6 = bsxfun(@and,tri2intpt6(:,v), tet2intpt6(:,tet_todo));
         common_intpts7 = bsxfun(@and,tri2intpt7(:,v), tet2intpt7(:,tet_todo));
         common_intpts8 = bsxfun(@and,tri2intpt8(:,v), tet2intpt8(:,tet_todo));
         common_intpts9 = bsxfun(@and,tri2intpt9(:,v), tet2intpt9(:,tet_todo));
         top_nodes_tet = bsxfun(@and,top_node2tet(:,tet_todo), rmdl.node2elem(:,v));
         bot_nodes_tet = bsxfun(@and,bot_node2tet(:,tet_todo), rmdl.node2elem(:,v));
      end
      tet_nodes     = bsxfun(@and,tet_node2tri(:,v), fmdl.node2elem(:,tet_todo));

      C = [C; v*ones(numel(tet_todo),1)];
      F = [F; tet_todo'];
      last_v = numel(V);
      V = [V; zeros(numel(tet_todo),1)]; % pre-allocate
      for t = 1:numel(tet_todo)
         pts = [  intpts1(common_intpts1(:,t),:);
                  intpts2(common_intpts2(:,t),:);
                  intpts3(common_intpts3(:,t),:);
                  fmdl.nodes(tet_nodes(:,t),:);];
         if ~isinf(opt.z_depth) 
            pts = [ pts;
                  intpts4(common_intpts4(:,t),:);
                  intpts5(common_intpts5(:,t),:);
                  intpts6(common_intpts6(:,t),:);
                  intpts7(common_intpts7(:,t),:);
                  intpts8(common_intpts8(:,t),:);
                  intpts9(common_intpts9(:,t),:);
                  top_nodes(top_nodes_tet(:,t),:);
                  bot_nodes(bot_nodes_tet(:,t),:)];
         end
         last_v = last_v + 1;
         
         if size(pts,1) < 4
%             debug_plot(fmdl,rmdl,v,tet_todo(t), bot, top, pts)
%             keyboard
            continue
         end
         if any(isnan(pts(:))), keyboard, end
%          debug_plot(fmdl,rmdl,v,tet_todo(t), bot, top, pts);
         try
            % move points to origin (helps for small elements at
            % large coordinates
            ctr = mean(pts);
            if any(isnan(ctr)), keyboard,end
            pts = bsxfun(@minus,pts,ctr);
            scale = max(abs(pts(:)));
            if scale == 0 %happens when there's only one point
               continue
            end
            % scale largest coordinate to 1 (helps with precision)
            pts = pts ./ scale;
            % force thorough search for initinal simplex and
            % supress precision warnings
            [K, V(last_v)] = convhulln(pts,{'Qt Pp Qs'});
            V(last_v) = max(V(last_v),0); % numerical issues may produce tiny negative volume
            V(last_v) = V(last_v) * scale^3; % undo scaling
         catch err
            ok = false;
            switch err.identifier
               case {'MATLAB:qhullmx:DegenerateData', 'MATLAB:qhullmx:UndefinedError'}
                  pts = bsxfun(@plus, pts .* scale, ctr);
                  u = uniquetol(pts,1e-14,'ByRows',true,'DataScale', 1);
                  ok = ok | size(u,1) < 4;
                  if ~ok
                     % test for colinearity in the xy plane
                     u12 = uniquetol(pts(:,1:2),1e-14,'ByRows',true,'DataScale',1);
                     cp = bsxfun(@minus, u12(2:end,1:2), u12(1,1:2));
                     l = sqrt(sum(cp.^2,2));
                     cp = bsxfun(@rdivide, cp, l);
                     % counteract colinear vectors in different directions
                     cp = abs(cp); 
                     un = uniquetol(cp,1e-12,'ByRows',true,'DataScale',1);
                     ok = ok | size(un,1) == 1; % co-linear points
                  end
                  if ~ok
                     % test if all points lie on the top or bottom caps
                     ok = ok | all(abs(pts(:,3) - top) < eps);
                     ok = ok | all(abs(pts(:,3) - bot) < eps);
                  end
            end
            if ~ok
               if DEBUG || eidors_debug('query','mk_tri2tet_c2f:convhulln')
                  debug_plot(fmdl,rmdl,v,tet_todo(t), bot, top, pts)
                  keyboard
               else
                  fprintf('\n');
                  eidors_msg(['convhulln has thrown an error. ' ...
                     'Enable eidors_debug on mk_tri2tet_c2f and re-run to see a debug plot'],0);
                  rethrow(err);
               end
            end
         end
      end
   end
   progress_msg(Inf);
   c2f = sparse(F,C,V,size(fmdl.elems,1),size(rmdl.elems,1));
   
   % add rtri contained in ftet
   if ~isinf(opt.z_depth)
      try rmdl = rmfield(rmdl,'coarse2fine'); end % messes with volume
      c2f = c2f + bsxfun(@times, sparse(tri_in_tet), opt.height * get_elem_volume(rmdl))';
   end
   % normalize to tet volume
   vol = get_elem_volume(fmdl);
   c2f = bsxfun(@rdivide,c2f,vol);
   
   % add tets contained in vox
   
   c2f = c2f + tet_in_tri;
   
end

function [intpts, tri2edge, tri2pts, edge2pts] = get_tet_intersection_pts(fmdl,rmdl,top,bot, epsilon)
   intpts = [];
   pt_idx = unique(rmdl.elems);
   pts = rmdl.nodes(pt_idx,1:2);
   
   bb_min = min(...
                min(fmdl.nodes(fmdl.faces(:,1),1:2),...
                    fmdl.nodes(fmdl.faces(:,2),1:2)),...
                fmdl.nodes(fmdl.faces(:,3),1:2));
   bb_max = max(...
                max(fmdl.nodes(fmdl.faces(:,1),1:2),...
                    fmdl.nodes(fmdl.faces(:,2),1:2)),...
                fmdl.nodes(fmdl.faces(:,3),1:2));
   todo = ~(  bsxfun(@gt,bb_min(:,1),pts(:,1)') ...
            | bsxfun(@gt,bb_min(:,2),pts(:,2)') ...
            | bsxfun(@lt,bb_max(:,1),pts(:,1)') ...
            | bsxfun(@lt,bb_max(:,2),pts(:,2)')); 
   [F,P] = find(todo);
   P = unique(P);
   in = false(size(fmdl.faces,1),size(pts,1));
   in(F,P) = point_in_triangle(pts(P,:),fmdl.faces(F,:),fmdl.nodes(:,1:2),epsilon)';
   
   [F,P] = find(in);

   % remove "vertical" faces
   vf = fmdl.normals(F,3) == 0;
   F(vf) = [];
   P(vf) = [];
   
   % project on faces
   % plane equation is ax+by+cz+d = 0, where d = -(ax0 + by0 + cz0)
   z = sum(fmdl.normals(F,:) .* fmdl.nodes(fmdl.faces(F,1),:),2);
%    z = repmat(d,1,length(P));
   for j = 1:2
      z = z - fmdl.normals(F,j) .* pts(P,j);
   end
   z = z ./ fmdl.normals(F,3);
   out = z>top | z < bot;
   F(out) = [];
   P(out) = [];
   z(out) = [];
   
   intpts = [pts(P,:) z];
   I = (1:size(intpts,1))';
   tri2edge = sparse(F,pt_idx(P),I,size(fmdl.faces,1),size(rmdl.nodes,1));
   tri2pts = sparse(F,I,ones(size(I,1),1),size(fmdl.faces,1),size(intpts,1));
   edge2pts = sparse(pt_idx(P),I,ones(size(I,1),1),size(rmdl.nodes,1),size(intpts,1));
   
end

function [intpts, FE2CE, FE2pts, CE2pts] = ...
   find_edge2edge_intersections_2d(FE, FN, CE, CN, top, bot)

   P1 = FN(FE(:,1),:);
   P2 = FN(FE(:,2),:);
   P3 = CN(CE(:,1),:);
   P4 = CN(CE(:,2),:);

   C_bb = zeros(size(P3,1),4);
   C_bb(:,[1 3]) = min(P3(:,1:2),P4(:,1:2));
   C_bb(:,[2 4]) = max(P3(:,1:2),P4(:,1:2));

   F_bb = zeros(size(P1,1),4);
   F_bb(:,[1 3]) = min(P1(:,1:2),P2(:,1:2));
   F_bb(:,[2 4]) = max(P1(:,1:2),P2(:,1:2));


   todo =   bsxfun(@gt, F_bb(:,1), C_bb(:,2)') ...
          | bsxfun(@lt, F_bb(:,2), C_bb(:,1)') ...
          | bsxfun(@gt, F_bb(:,3), C_bb(:,4)') ...
          | bsxfun(@lt, F_bb(:,4), C_bb(:,3)');
   todo = ~todo;

   [T, V] = find(todo);
   
   S1 = P2(T,:) - P1(T,:);
   S2 = P4(V,:) - P3(V,:);
   
   denom =    S2(:,2) .* S1(:,1) - S2(:,1) .* S1(:,2);

   P13 = P1(T,1:2) - P3(V,1:2);

   num_a =    S2(:,1) .* P13(:,2) ...
            - S2(:,2) .* P13(:,1);
   num_b =    S1(:,1) .* P13(:,2) ...
            - S1(:,2) .* P13(:,1);
   
   mua = num_a ./ denom;
   mub = num_b ./ denom;
   
   IN = mua>0 & mua<1 & mub>0 & mub<1;
   T = T(IN);
   V = V(IN);
   intpts = P1(T,:) + bsxfun(@times, mua(IN), S1(IN,:));
   in = ~(intpts(:,3) > top | intpts(:,3) < bot);
   intpts = intpts(in,:);
   I = (1:size(intpts,1))';
   T = T(in);
   V = V(in);
   FE2CE = sparse(size(P1,1),size(P3,1));
   FE2CE(sub2ind(size(FE2CE),T,V)) = I;
   FE2pts = sparse(T,I,ones(size(I)),size(P1,1),size(I,1));
   CE2pts  = sparse(V,I,ones(size(I)),size(P3,1),size(I,1));
   
end

function [top_node2tet, top_nodes, nodes_above, nodes_below] = ...
   get_cap_nodes_in_tets(fmdl,rmdl,top,epsilon)

   top_nodes = rmdl.nodes;
   top_nodes(:,3) = top;
   nodes_above = fmdl.nodes(:,3) >= top;
   nodes_below = fmdl.nodes(:,3) <= top;
   
   % nodes in tets
   elidx = any(nodes_above(fmdl.elems),2) & any(nodes_below(fmdl.elems),2);
   top_node2tet = sparse(size(rmdl.nodes,1),size(fmdl.elems,1));
   if any(elidx)
      mdl = struct;
      mdl.nodes = fmdl.nodes;
      mdl.elems = fmdl.elems(elidx,:);
      top_node2tet(:,elidx) = point_in_tet(mdl,top_nodes,epsilon);
   end
end

%-------------------------------------------------------------------------%
% Intersections between tet edges and tri faces
function [intpts, edge2tri, edge2pts, tri2pts] = ...
         get_cap_tet_edge2tri_face_intersections(fmdl,rmdl,top,...
                                          nodes_above,nodes_below, epsilon)
   
   intpts = [];
   edge2tri = sparse(size(fmdl.edges,1),size(rmdl.elems,1));
   edge2pts = sparse(size(fmdl.edges,1),0);
   tri2pts  = sparse(size(rmdl.elems,1),0);
   

   % tet_edge2tri_face
   edgeidx = any(nodes_above(fmdl.edges),2) & any(nodes_below(fmdl.edges),2);
   % discard nodes on the plane
   edgeidx(edgeidx) = fmdl.nodes(fmdl.edges(edgeidx,1),3) ~= ...
                      fmdl.nodes(fmdl.edges(edgeidx,2),3);
   
   v = fmdl.nodes(fmdl.edges(edgeidx,2),:) - ...
       fmdl.nodes(fmdl.edges(edgeidx,1),:);
   u = (top - fmdl.nodes(fmdl.edges(edgeidx,1),3))  ./ v(:,3);
   pts = fmdl.nodes(fmdl.edges(edgeidx,1),:) + bsxfun(@times,u, v);
   t = point_in_triangle(pts(:,1:2),rmdl.elems,rmdl.nodes(:,1:2),-epsilon);
   % remove any edges that cross multiple faces -- they will be cought by 
   % edge2edge intersections elsewhere
   t(sum(t,2)>1,:) = false;
   [E,T] = find(t);
   
   if any(t(:))
      intpts = pts(E,:);
      N = size(intpts,1);
      I = (1:N)';
      emap = find(edgeidx);
      E = emap(E);
      edge2tri = sparse(E,T,I,size(fmdl.edges,1), size(rmdl.elems,1));
      edge2pts = sparse(E,I,ones(size(I)),size(fmdl.edges,1), N);
      tri2pts  = sparse(T,I,ones(size(I)),size(rmdl.elems,1), N);
   end
   
end
%-------------------------------------------------------------------------%
% Intersections between tet faces and tri edges
function   [intpts, tri2edge,tri2intpt,edge2intpt] = ...
                           get_cap_tet_face2tri_edge_intersections( ...
                            fmdl,rmdl,top,nodes_above,nodes_below, epsilon)
   
   USE_POINT_IN_TRIANGLE = 0;
   
   intpts = [];
   tri2edge = sparse(size(fmdl.faces,1),size(rmdl.edges,1));
   tri2intpt = sparse(size(fmdl.faces,1),0);
   edge2intpt  = sparse(size(rmdl.edges,1),0);

                              
   %faces that cross the cap
   faceidx = any(nodes_above(fmdl.faces),2) & any(nodes_below(fmdl.faces),2);
 
   faces = fmdl.faces(faceidx,:);
   normals = fmdl.normals(faceidx,:);
   
   N_edges = size(rmdl.edges,1);
   N_faces = size(faces,1);
   
  
   face_bb = zeros(N_faces,6);
   face_bb(:,1) = min(reshape(fmdl.nodes(faces,1),N_faces,3),[],2);
   face_bb(:,2) = max(reshape(fmdl.nodes(faces,1),N_faces,3),[],2);
   face_bb(:,3) = min(reshape(fmdl.nodes(faces,2),N_faces,3),[],2);
   face_bb(:,4) = max(reshape(fmdl.nodes(faces,2),N_faces,3),[],2);
   
   edge_bb = zeros(N_edges,6);
   edge_bb(:,1) = min(reshape(rmdl.nodes(rmdl.edges,1),N_edges,2),[],2);
   edge_bb(:,2) = max(reshape(rmdl.nodes(rmdl.edges,1),N_edges,2),[],2);
   edge_bb(:,3) = min(reshape(rmdl.nodes(rmdl.edges,2),N_edges,2),[],2);
   edge_bb(:,4) = max(reshape(rmdl.nodes(rmdl.edges,2),N_edges,2),[],2);
   
   todo =   bsxfun(@gt, face_bb(:,1), edge_bb(:,2)') ...
          | bsxfun(@lt, face_bb(:,2), edge_bb(:,1)') ...
          | bsxfun(@gt, face_bb(:,3), edge_bb(:,4)') ...
          | bsxfun(@lt, face_bb(:,4), edge_bb(:,3)');
   todo = ~todo;
   e_todo = any(todo,1);
   f_todo = any(todo,2);
   faceidx(faceidx) = f_todo;
   faces = faces(f_todo,:);
   normals = normals(f_todo,:);
   [F,E] = find(todo);
   emap = zeros(size(e_todo,2),1);
   emap(e_todo) = 1:nnz(e_todo);
   E = emap(E);
   fmap = zeros(size(f_todo,1),1);
   fmap(f_todo) = 1:nnz(f_todo);
   F = fmap(F);
   
   P1 = rmdl.nodes(rmdl.edges(e_todo,1),:);
   P12 = P1 - rmdl.nodes(rmdl.edges(e_todo,2),:);
   P1(:,3) = top;
   P12(:,3) = 0;

   d = sum(fmdl.normals(faceidx,:) .* fmdl.nodes(fmdl.faces(faceidx,1),:),2);
   
   if ~USE_POINT_IN_TRIANGLE
      % for point_in_triangle
      nodes1 = fmdl.nodes(faces(:,1),:);
      v0 = fmdl.nodes(faces(:,3),:) - nodes1;
      v1 = fmdl.nodes(faces(:,2),:) - nodes1;
      dot00 = dot(v0, v0, 2);
      dot01 = dot(v0, v1, 2);
      % dot02 = dot(v0, v2, 2);
      dot11 = dot(v1, v1, 2);
      % dot12 = dot(v1, v2, 2);
      invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
   end
   
   % find points of intersection between edges and face planes
   num = -d(F) + sum(normals(F,:).*P1(E,:),2);
   den = sum(normals(F,:) .* P12(E,:),2);
   u = num ./ den;
   % den == 0 => normal perpendicular to line
   idx = u >= 0 & u <= 1 & abs(den) >= eps;
   
   if any(idx)      
      F = F(idx);
      E = E(idx);
      ipts = P1(E,:) - bsxfun(@times, u(idx), P12(E,:));
      if USE_POINT_IN_TRIANGLE
         t = point_in_triangle(ipts,faces(F,:),fmdl.nodes,epsilon,'match');
      else
         v2 = ipts - fmdl.nodes(faces(F,1),:);
         dot02 = dot(v0(F,:),v2,2);
         dot12 = dot(v1(F,:),v2,2);
         % barycentric coordinates
         dot01 = dot01(F);
         invDenom = invDenom(F);
         u = (dot11(F) .* dot02 - dot01 .* dot12) .* invDenom;
         v = (dot00(F) .* dot12 - dot01 .* dot02) .* invDenom;
         t = u >= -epsilon & v >= -epsilon & (u+v-epsilon) <= 1;
      end
      
      if any(t)
         N = nnz(t);
         idv = (1:N)';
         intpts = ipts(t,:);
         I = idv;
         idx = find(faceidx); idx = idx(F);
         F = idx(t);
         eimap = find(emap); 
         E = eimap(E(t));
         
         tri2edge = sparse(F,E,I,size(fmdl.faces,1),size(rmdl.edges,1));
         tri2intpt = sparse(F,I,ones(size(I)),size(fmdl.faces,1),size(I,1));
         edge2intpt  = sparse(E,I,ones(size(I)),size(rmdl.edges,1),size(I,1));
      end
   end
end

%-------------------------------------------------------------------------%
% Center scale models
function [fmdl,rmdl, opt] = center_scale_models(fmdl,rmdl, opt)
   ctr = mean([min(rmdl.nodes);max(rmdl.nodes)]);
   rmdl.nodes = bsxfun(@minus,rmdl.nodes,ctr);
   fmdl.nodes = bsxfun(@minus,fmdl.nodes,ctr);
   opt.top = opt.top - ctr(3);
   opt.bot = opt.bot - ctr(3);
   if isfield(opt,'do_not_scale') && opt.do_not_scale
      return
   end
   maxnode = min( max(abs(rmdl.nodes(:))), max(abs(fmdl.nodes(:))));
   scale = 1/maxnode;
   rmdl.nodes = scale*rmdl.nodes;
   fmdl.nodes = scale*fmdl.nodes;
   opt.top    = scale*opt.top;
   opt.bot    = scale*opt.bot;
   opt.height = scale*opt.height;
   opt.z_depth= scale*opt.z_depth;
   eidors_msg('@@ models scaled by %g', scale,2);
end

%-------------------------------------------------------------------------%
% Remove obviously non-overlapping parts of the models
function [fmdl,rmdl,fmdl_idx,rmdl_idx] = crop_models(fmdl,rmdl, opt)
   f_min = min(fmdl.nodes);
   f_max = max(fmdl.nodes);
   r_min = min(rmdl.nodes);
   r_max = max(rmdl.nodes);
   
   if isfield(opt, 'z_depth') && ~isinf(opt.z_depth)
      lvl = mean(rmdl.nodes(:,3));
      r_max(3) = lvl + opt.z_depth;
      r_min(3) = lvl - opt.z_depth;
   else
      r_max(3) = f_max(3);
      r_min(3) = f_min(3);
   end
   
   % nodes outside the bounding box of the other model
   f_gt  = bsxfun(@gt, fmdl.nodes, r_max);
   f_lt  = bsxfun(@lt, fmdl.nodes, r_min);
   r_gt  = bsxfun(@gt, rmdl.nodes, f_max);
   r_lt  = bsxfun(@lt, rmdl.nodes, f_min);
   
   % elems outside the bounding box of the other model
   re_gt = any(reshape(all(reshape(r_gt(rmdl.elems',:),3,[])),[],3),2);
   re_lt = any(reshape(all(reshape(r_lt(rmdl.elems',:),3,[])),[],3),2);
   fe_gt = any(reshape(all(reshape(f_gt(fmdl.elems',:),4,[])),[],3),2);
   fe_lt = any(reshape(all(reshape(f_lt(fmdl.elems',:),4,[])),[],3),2);
   
   % elems to keep
   rmdl_idx = ~(re_gt | re_lt);
   fmdl_idx = ~(fe_gt | fe_lt);
   
   % remove non-overlapping elems
   rmdl.elems = rmdl.elems(rmdl_idx,:);
   fmdl.elems = fmdl.elems(fmdl_idx,:);
   
   % remove unused nodes
   [r_used_nodes,jnk,r_n] = unique(rmdl.elems(:));
   [f_used_nodes,jnk,f_n] = unique(fmdl.elems(:));
   
   r_idx = 1:numel(r_used_nodes);
   f_idx = 1:numel(f_used_nodes);
   
   rmdl.elems = reshape(r_idx(r_n),[],3);
   fmdl.elems = reshape(f_idx(f_n),[],4);
   
   rmdl.nodes = rmdl.nodes(r_used_nodes,:);
   fmdl.nodes = fmdl.nodes(f_used_nodes,:);
   
   % for the benefit of any (debug) plots later on
   try, rmdl = rmfield(rmdl,'boundary'); end
   try, fmdl = rmfield(fmdl,'boundary'); end
end

%-------------------------------------------------------------------------%
% Prepare model
function fmdl = prepare_tet_mdl(fmdl)
   fmopt.elem2edge = true;
   fmopt.edge2elem = true;
   fmopt.face2elem = true;
   fmopt.node2elem = true;
   fmopt.normals   = true;
   ll = eidors_msg('log_level',1);
   fmopt.linear_reorder = false; % this is slow and not needed
   fmdl = fix_model(fmdl,fmopt);
   eidors_msg('log_level',ll);
   fmdl.node2elem = logical(fmdl.node2elem);
   nElem = size(fmdl.elems,1);
   nFace = size(fmdl.faces,1);
   fmdl.elem2face = sparse(repmat((1:nElem)',1,4),double(fmdl.elem2face),true,nElem,nFace);
end

%-------------------------------------------------------------------------%
% Prepare model
function fmdl = prepare_tri_mdl(fmdl)
   fmopt.elem2edge = true;
   fmopt.edge2elem = true;
%    fmopt.face2elem = true;
   fmopt.node2elem = true;
%    fmopt.normals   = true;
   fmopt.linear_reorder = false; % this is slow and not needed
   ll = eidors_msg('log_level',1);
   fmdl = fix_model(fmdl,fmopt);
   eidors_msg('log_level',ll);
   fmdl.node2elem = logical(fmdl.node2elem);

end

%-------------------------------------------------------------------------%
% Debug plot
function debug_plot(fmdl,rmdl,v,t, bot, top, pts)
   tet.nodes = fmdl.nodes;
   tri.nodes = repmat(rmdl.nodes(rmdl.elems(v,:),:),2,1);
   tri.nodes(:,3) = [repmat(bot,3,1); repmat(top,3,1)];
   tri.elems = [ 1 2 5 4
                 2 3 6 5
                 3 1 4 6];
   tri.boundary = tri.elems;
   tet.type = 'fwd_model';
   tri.type = 'fwd_model';
   tet.elems = fmdl.elems(t,:);
   clf
   show_fem(tri)
   hold on
   h = show_fem(tet);
   set(h,'EdgeColor','b')
%    pts = bsxfun(@plus,pts*scale,ctr);
   plot3(pts(:,1),pts(:,2),pts(:,3),'o');
   hold off
   axis auto
end

%-------------------------------------------------------------------------%
% Parse option struct
function opt = parse_opts(fmdl,rmdl, opt)

   lvl = mean(rmdl.nodes(:,3));
   
   if ~isfield(opt, 'z_depth')
      opt.z_depth = Inf;
   end
   if isinf(opt.z_depth)
      opt.top = max(fmdl.nodes(:,3));
      opt.bot = min(fmdl.nodes(:,3));
      opt.height = opt.top - opt.bot;
   else
      opt.top = lvl + opt.z_depth;
      opt.bot = lvl - opt.z_depth;
      opt.height = 2 * opt.z_depth;
   end
   if ~isfield(opt, 'tol_node2tri');
      opt.tol_node2tri = eps; % * max(rmdl_rng,fmdl_rng)^3;
   end
   if ~isfield(opt, 'tol_node2tet');
      opt.tol_node2tet = eps; % * max(rmdl_rng,fmdl_rng)^3;
   end
   if ~isfield(opt, 'tol_edge2edge')
      opt.tol_edge2edge = 6*eps(min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))));
   end
   if ~isfield(opt, 'tol_edge2tri')
      opt.tol_edge2tri = eps; %1e-10
   end
   %     if ~isfield(opt, 'save_memory')
   %        opt.save_memory = 0;
   %     end
   eidors_msg('@@ node2tet  tolerance = %g', opt.tol_node2tet,2);
   eidors_msg('@@ edge2edge tolerance = %g', opt.tol_edge2edge,2);
   eidors_msg('@@ edge2tri  tolerance = %g', opt.tol_edge2tri,2);
 end
 
 
function do_unit_test
   do_small_test
   do_inf_test
   do_realistic_test
   do_centre_slice; % failing code from tutorial
end

function do_inf_test
   jnk = mk_common_model('a2c',16);
   tri = jnk.fwd_model;
   tri.nodes = 1.1*tri.nodes;
   jnk = mk_common_model('c3cr',16);
   tet = jnk.fwd_model;

   c2f_a = mk_tri2tet_c2f(tet,tri);
   c2f_n = mk_approx_c2f(tet,tri);
   
   prob = c2f_n ~= 0 & c2f_a == 0;
   unit_test_cmp('Assert no missing intersections', nnz(prob),0);
   unit_test_cmp('mk_tri2tet_c2f v approximate', c2f_n,c2f_a, .5);
end


function do_small_test
   jnk = mk_common_model('a2c',16);
   tri = jnk.fwd_model;
   tri.nodes = 1.1*tri.nodes;
   jnk = mk_common_model('c3cr',16);
   tet = jnk.fwd_model;
%    tet.elems = tet.elems(8081,:);
%    tri.elems = tri.elems(1,:);
   opt.z_depth = 0.5;
   c2f = mk_tri2tet_c2f(tet,tri,opt);
 
   clf
   subplot(131)
   show_fem(tet);
   hold on
   show_fem(tri);
   view(2)
   subplot(132)
   img = mk_image(tet,1);
   img.elem_data = sum(c2f,2);
   img.calc_colours.clim = 1;
   img.calc_colours.ref_level = 0;
   show_fem(img);
   subplot(133)
   img = mk_image(tri,1);
   img.elem_data = sum(bsxfun(@times,c2f,get_elem_volume(tet)),1) ./ get_elem_volume(tri)';
   img.calc_colours.clim = 1;
   img.calc_colours.ref_level = 0;
   show_fem(img)
   
   unit_test_cmp('Check C2F size  ', size(c2f),[length(tet.elems), length(tri.elems)]);
   unit_test_cmp('Check C2F max==1', max(c2f(:)), 1);
   unit_test_cmp('Check C2F min==0', min(c2f(:)), 0);
   
   f2c = bsxfun(@rdivide, bsxfun(@times, c2f, get_elem_volume(tet))', get_elem_volume(tri));
   unit_test_cmp('Check F2C max==1', max(sum(f2c,2)), 1, 1e-14);
   unit_test_cmp('Check F2C min==0', min(f2c(:)), 0);
end
 

function do_realistic_test
      
   fmdl= ng_mk_cyl_models([2,2,.1],[16,1],[.1,0,.025]);
   xvec = -2:.5:2; %[-1.5 -.5:.2:.5 1.5];
   yvec = -2:.5:2; %[-1.6 -1:.2:1 1.6];
   zvec = [.5 1.5];
   grid2d = mk_grid_model([],xvec,yvec);
   grid3d = mk_grid_model([],xvec,yvec,zvec);
   eidors_cache clear
   tic
   opt.save_memory = 0;
   c2f_a = mk_grid_c2f(fmdl, grid3d,opt);
   t = toc;
   fprintf('Voxel:\tt=%f s\n',t);

   opt.z_depth = .5;
   grid2d.nodes(:,3) = 1;
   eidors_cache clear
   tic
%    fmdl.elems = fmdl.elems(2764,:);
%    grid2d.elems = grid2d.elems(4,:);   
   c2f_b = mk_tri2tet_c2f(fmdl, grid2d, opt);
   t = toc;
   fprintf('tri2tet:\tt=%f \n',t);
   c2f_c = c2f_b * grid2d.coarse2fine;
   
   eidors_cache clear

   grid2d.mk_coarse_fine_mapping.z_depth = .5;
   grid2d.mk_coarse_fine_mapping.f2c_offset(3) = 1;
   grid2d = rmfield(grid2d,'coarse2fine');
   tic
   c2f_n = mk_approx_c2f(fmdl,grid2d);
   t = toc;
   fprintf('Approximate: t=%f s\n',t);
   
   unit_test_cmp('mk_tri2tet_c2f v mk_grid_c2f', c2f_c,c2f_a, 1e-5);
   unit_test_cmp('mk_tri2tet_c2f v approximate', c2f_n,c2f_b, .5);
   prob = c2f_n ~= 0 & c2f_b == 0;
   unit_test_cmp('Assert no missing intersections', nnz(prob),0);
%    keyboard
%    subplot(132)
%    imagesc(c2f_c, [0 1]);
%    subplot(133)
%    imagesc(c2f_a, [0 1]);
%    figure
%    subplot(131)
% show_fem(mk_image(fmdl,sum(c2f_a,2)));
% img_a = mk_image(fmdl,sum(c2f_a,2));
% img_a.calc_colours.ref_level = 0;
% show_fem(img_a)
% subplot(132)
% img_c = mk_image(fmdl,sum(c2f_c,2));
% img_c.calc_colours.ref_level = 0;
% show_fem(img_c)
% subplot(133)
% img_n = mk_image(fmdl,sum(c2f_n,2));
% img_n.calc_colours.ref_level = 0;
% show_fem(img_n)
%    keyboard
end
 
 
function do_centre_slice; % failing code from tutorial
   imdl = mk_common_model('b3cr',[16,2]);
   f_mdl= imdl.fwd_model;
   imdl2d= mk_common_model('b2c2',16);
   c_mdl= imdl2d.fwd_model;
   c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
end
