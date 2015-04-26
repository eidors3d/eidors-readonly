function [c2f] = mk_tet_c2f(fmdl, rmdl, opt)
%MK_TET_C2F - calculate a coarse2fine mapping for two tet-based models.
% C2F = MK_TET_C2F(FMDL,RMDL) returns in C2F the fraction of volume of
% each element of the fine model FMDL contained in each element of
% the coarse model RMDL.
% Uses CONVHULLN to calculate the volume defined by a set of intersection
% points between individual tet and vox elements.
%
% C2F = MK_TET_C2F(FMDL,RMDL,OPT) allows specifying options.
% 
% Inputs:
%   FMDL - a (fine) EIDORS (tet-based) forward model
%   RMDL - a (course) EIDORS (tet-based) forward model
%   OPT  - an option structure with the following fields and defaults:
%      .do_not_scale  - set to true to prevent scaling the models to unit
%                       cube before any calculations, including thresholds.
%                       Default: false
%      .tol_node2tet  - tolerance for determinant <= 0 in testing for
%                       points inside tets. Default: eps
%      .tol_edge2edge - maximum distance between "intersecting" edges
%                       Default: 6*sqrt(3)*eps(a), where a is
%                       min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))
%      .tol_edge2tri  - minimum value of a barycentric coordinate to 
%                       decide a point is lying inside a triangle and not
%                       on its edge. Default: eps
%
% NOTE that for grid-based models, such as returned by MK_GRID_MODEL or
% MK_VOXEL_VOLUME, MK_GRID_C2F is much faster.
%
% Set eidors_msg 'log level' < 2 to supress output to command line.
%
% Examples:
%     fmdl = ng_mk_cyl_models([2,2,.2],[],[]);
%     rmdl = ng_mk_cyl_models([2,2],[],[]);
%     c2f  = mk_tet_c2f(fmdl,rmdl);
%     h = show_fem(fmdl); set(h,'LineWidth',0.1)
%     hold on
%     h = show_fem(rmdl); set(h,'EdgeColor','b','LineWidth',2);
%     hold off
%
% See also MK_GRID_C2F, FIND_EDGE2EDGE_INTERSECTIONS, CONVHULLN
%          MK_COARSE_FINE_MAPPING, POINT_IN_TRIANGLE, EIDORS_MSG


% (C) 2015 Bartlomiej Grychtol - all rights reserved by Swisstom AG
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'), do_unit_test; return; end
if nargin < 3
   opt = struct();
end

f_elems = size(fmdl.elems,1);
r_elems = size(rmdl.elems,1);

c2f = sparse(f_elems,r_elems);
[fmdl,rmdl,fmdl_idx,rmdl_idx] = crop_models(fmdl,rmdl);

if ~(any(fmdl_idx) && any(rmdl_idx))
   eidors_msg('@@: models do not overlap, returning all-zeros');
   return
end

[fmdl,rmdl] = center_scale_models(fmdl,rmdl, opt);

opt = parse_opts(fmdl,rmdl, opt);


copt.fstr = 'mk_tet_c2f';

c2f(fmdl_idx,rmdl_idx) = eidors_cache(@do_mk_tet_c2f,{fmdl,rmdl,opt},copt);


function c2f = do_mk_tet_c2f(fmdl,rmdl,opt)
   DEBUG = eidors_debug('query','mk_tet_c2f');
   
   c2f = sparse(0,0);
   progress_msg('Prepare fine model...');
   fmdl = prepare_tet_mdl(fmdl);
   progress_msg(Inf);
   
   progress_msg('Prepare course model...');
   rmdl = prepare_tet_mdl(rmdl);
   progress_msg(Inf);
   
   progress_msg('Find c_edge2f_face intersections...')
   [intpts1, fface2redge, fface2intpt1, redge2intpt1] = ...
      edge2face_intersections(fmdl,rmdl,opt);
   progress_msg(sprintf('Found %d', size(intpts1,1)), Inf);

   progress_msg('Find f_edge2c_face intersections...')
   [intpts2, rface2fedge, rface2intpt2, fedge2intpt2] = ...
      edge2face_intersections(rmdl,fmdl,opt);
   progress_msg(sprintf('Found %d', size(intpts2,1)), Inf);

   pmopt.final_msg = 'none';
   progress_msg('Find edge2edge intersections...',-1,pmopt)
   [intpts3, fedge2redge, fedge2intpt3, redge2intpt3] = ...
      find_edge2edge_intersections(fmdl.edges, fmdl.nodes, ...
                                   rmdl.edges, rmdl.nodes, ...
                                   opt.tol_edge2edge);
   progress_msg(sprintf('Found %d',size(intpts3,1)),Inf);

   progress_msg('Find c_nodes in f_tets...');
   rnode2ftet = get_nodes_in_tets(fmdl,rmdl.nodes, opt);
   progress_msg(sprintf('Found %d', nnz(rnode2ftet)), Inf);
   
   
   progress_msg('Find c_elems in f_elems...')
   rtet_in_ftet = (double(rmdl.node2elem') * rnode2ftet) == 4;
   progress_msg(sprintf('Found %d',nnz(rtet_in_ftet)), Inf);
   
   progress_msg('Find f_nodes in c_tets...');
   fnode2rtet = get_nodes_in_tets(rmdl,fmdl.nodes, opt);
   progress_msg(sprintf('Found %d', nnz(fnode2rtet)), Inf);

   progress_msg('Find f_elems in c_elems...')
   ftet_in_rtet = (double(fmdl.node2elem') * fnode2rtet) == 4;
   progress_msg(sprintf('Found %d',nnz(ftet_in_rtet)), Inf);
   
   progress_msg('Find total intersections...');
   e2e = double(rmdl.edge2elem');
   rtet2ftet =  double(rmdl.elem2face) * (rface2fedge>0) * fmdl.edge2elem ...
                 | e2e * (fface2redge>0)' * fmdl.elem2face' ...
                 | e2e * fedge2redge' * fmdl.edge2elem;
   % exclude inclusion (dealt with separately)
   rtet2ftet = rtet2ftet & ~rtet_in_ftet & ~ftet_in_rtet'; 
   progress_msg(sprintf('Found %d',nnz(rtet2ftet)), Inf);

   
   progress_msg('Calculate intersection volumes...');
   % sparse logical multiplication doesn't exist
   rtet2intpt1 = logical(rmdl.edge2elem'*redge2intpt1)';
   ftet2intpt1 = logical(fmdl.elem2face *fface2intpt1)';
   
   rtet2intpt2 = logical(rmdl.elem2face * rface2intpt2)';
   ftet2intpt2 = logical(fmdl.edge2elem'* fedge2intpt2)';
   
   ftet2intpt3 = logical(fmdl.edge2elem'* fedge2intpt3)';
   rtet2intpt3 = logical(rmdl.edge2elem'* redge2intpt3)';
    
   rtet_todo = find(sum(rtet2ftet,2)>0);
   C = []; F = []; V = [];
   
   id = 0; N = length(rtet_todo);
   mint = ceil(N/100);
   for v = rtet_todo'
      id = id+1;
      if mod(id,mint)==0, progress_msg(id/N); end
      tet_todo = find(rtet2ftet(v,:));
      common_intpts1 = bsxfun(@and,rtet2intpt1(:,v), ftet2intpt1(:,tet_todo));
      common_intpts2 = bsxfun(@and,rtet2intpt2(:,v), ftet2intpt2(:,tet_todo));
      common_intpts3 = bsxfun(@and,rtet2intpt3(:,v), ftet2intpt3(:,tet_todo));
      f_nodes     = bsxfun(@and,fnode2rtet(:,v), fmdl.node2elem(:,tet_todo));
      r_nodes     = bsxfun(@and,rnode2ftet(:,tet_todo), rmdl.node2elem(:,v));
      C = [C; v*ones(numel(tet_todo),1)];
      F = [F; tet_todo'];
      last_v = numel(V);
      V = [V; zeros(numel(tet_todo),1)]; % pre-allocate
      
      for t = 1:numel(tet_todo)
         pts = [ intpts1(common_intpts1(:,t),:);
            intpts2(common_intpts2(:,t),:);
            intpts3(common_intpts3(:,t),:);
            fmdl.nodes(f_nodes(:,t),:);
            rmdl.nodes(r_nodes(:,t),:)];
         last_v = last_v + 1;
         try
            % move points to origin (helps for small elements at
            % large coordinates
            ctr = mean(pts);
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
            V(last_v) = V(last_v) * scale^3; % undo scaling
         catch err
            ok = false;
            switch err.identifier
               case {'MATLAB:qhullmx:DegenerateData', 'MATLAB:qhullmx:UndefinedError'}
                  if size(pts,1) > 3
                     u = uniquetol(pts*scale,eps,'ByRows',true,'DataScale', 1);
                     ok = ok | size(u,1) < 4;
                  end
            end
            if ~ok
               if DEBUG || eidors_debug('query','mk_tet_c2f:convhulln');
                  tet.nodes = fmdl.nodes;
                  vox.nodes = rmdl.nodes;
                  tet.type = 'fwd_model';
                  vox.type = 'fwd_model';
                  vox.elems = rmdl.faces(logical(rmdl.elem2face(v,:)),:);
                  vox.boundary = vox.elems;
                  tet.elems = fmdl.elems(tet_todo(t),:);
                  clf
                  show_fem(vox)
                  hold on
                  h = show_fem(tet);
                  set(h,'EdgeColor','b')
                  pts = bsxfun(@plus,pts*scale,ctr);
                  plot3(pts(:,1),pts(:,2),pts(:,3),'o');
                  hold off
                  axis auto
                  keyboard
               else
                  fprintf('\n');
                  eidors_msg(['convhulln has thrown an error. ' ...
                     'Enable eidors_debug on mk_tet_c2f and re-run to see a debug plot'],0);
                  rethrow(err);
               end
            end
         end
      end
   end
   progress_msg(Inf);
    
    c2f = sparse(F,C,V,size(fmdl.elems,1),size(rmdl.elems,1));
    
    % add rtet contained in ftet
    try rmdl = rmfield(rmdl,'coarse2fine'); end % messes with volume
    c2f = c2f + bsxfun(@times, sparse(rtet_in_ftet), get_elem_volume(rmdl))';
    
    % normalize to tet volume
    vol = get_elem_volume(fmdl);
    c2f = bsxfun(@rdivide,c2f,vol);

    % add tets contained in vox

    c2f = c2f + ftet_in_rtet;
%-------------------------------------------------------------------------%
% Calculate intersection points between faces and edges    
function [intpts, tri2edge, tri2intpt, edge2intpt] = edge2face_intersections(fmdl,rmdl,opt)
   N_edges = size(rmdl.edges,1);
   N_faces = size(fmdl.faces,1);
   
   face_bb = zeros(N_faces,6);
   face_bb(:,1) = min(reshape(fmdl.nodes(fmdl.faces,1),N_faces,3),[],2);
   face_bb(:,2) = max(reshape(fmdl.nodes(fmdl.faces,1),N_faces,3),[],2);
   face_bb(:,3) = min(reshape(fmdl.nodes(fmdl.faces,2),N_faces,3),[],2);
   face_bb(:,4) = max(reshape(fmdl.nodes(fmdl.faces,2),N_faces,3),[],2);
   face_bb(:,5) = min(reshape(fmdl.nodes(fmdl.faces,3),N_faces,3),[],2);
   face_bb(:,6) = max(reshape(fmdl.nodes(fmdl.faces,3),N_faces,3),[],2);
   
   edge_bb = zeros(N_edges,6);
   edge_bb(:,1) = min(reshape(rmdl.nodes(rmdl.edges,1),N_edges,2),[],2);
   edge_bb(:,2) = max(reshape(rmdl.nodes(rmdl.edges,1),N_edges,2),[],2);
   edge_bb(:,3) = min(reshape(rmdl.nodes(rmdl.edges,2),N_edges,2),[],2);
   edge_bb(:,4) = max(reshape(rmdl.nodes(rmdl.edges,2),N_edges,2),[],2);
   edge_bb(:,5) = min(reshape(rmdl.nodes(rmdl.edges,3),N_edges,2),[],2);
   edge_bb(:,6) = max(reshape(rmdl.nodes(rmdl.edges,3),N_edges,2),[],2);
   
   allocsz = max(N_edges,N_faces);
   N_alloc = allocsz;
   
   intpts = zeros(N_edges,3);
   T = zeros(N_edges,1);
   E = zeros(N_edges,1);
   I = zeros(N_edges,1);
  
   P1 = rmdl.nodes(rmdl.edges(:,1),:);
   P12 = P1 - rmdl.nodes(rmdl.edges(:,2),:);

   
   d = sum(fmdl.normals .* fmdl.nodes(fmdl.faces(:,1),:),2);
      
   mint = ceil(N_edges/100);
   
   % for point_in_triangle
   v0 = fmdl.nodes(fmdl.faces(:,3),:) - fmdl.nodes(fmdl.faces(:,1),:);
   v1 = fmdl.nodes(fmdl.faces(:,2),:) - fmdl.nodes(fmdl.faces(:,1),:);
   dot00 = dot(v0, v0, 2);
   dot01 = dot(v0, v1, 2);
   % dot02 = dot(v0, v2, 2);
   dot11 = dot(v1, v1, 2);
   % dot12 = dot(v1, v2, 2);
   invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
   
   epsilon = opt.tol_edge2tri;
   
   excl =   bsxfun(@gt, face_bb(:,1), edge_bb(:,2)') ...
          | bsxfun(@lt, face_bb(:,2), edge_bb(:,1)') ...
          | bsxfun(@gt, face_bb(:,3), edge_bb(:,4)') ...
          | bsxfun(@lt, face_bb(:,4), edge_bb(:,3)') ...
          | bsxfun(@gt, face_bb(:,5), edge_bb(:,6)') ...
          | bsxfun(@lt, face_bb(:,6), edge_bb(:,5)');
   excl = ~excl;
   N_pts = 0;
   for i = 1:N_edges
      if mod(i,mint)==0, progress_msg(i/N_edges); end
     
      fidx = excl(:,i);
      if ~any(fidx), continue, end;
      
      num = -d(fidx) + sum(bsxfun(@times,fmdl.normals(fidx,:),P1(i,:)),2);
      
      den = sum(bsxfun(@times,fmdl.normals(fidx,:),P12(i,:)),2);
      
      u = num ./ den;
      
      idx = u >= 0 & u <= 1;
      
      % calculate the intersection points
      if any(idx)
         id = find(idx);
         ipts = bsxfun(@minus, P1(i,:), bsxfun(@times, u(id), P12(i,:)));
         
         if 1
            fcs = find(fidx);
            fid = fcs(id);
            % point in triangle test
            v2 = bsxfun(@minus,ipts,fmdl.nodes(fmdl.faces(fid,1),:));
            dot02 = dot(v0(fid,:),v2,2);
            dot12 = dot(v1(fid,:),v2,2);
            % barycentric coordinates
            u = (dot11(fid) .* dot02 - dot01(fid) .* dot12) .* invDenom(fid);
            v = (dot00(fid) .* dot12 - dot01(fid) .* dot02) .* invDenom(fid);
            t = u >= -epsilon & v >= -epsilon & (u+v-epsilon) <= 1; 
         else
            t = point_in_triangle(ipts,fmdl.faces(id,:),fmdl.nodes,epsilon,'match');
         end
         if any(t)
            N = nnz(t);
            if N_pts+N > N_alloc
               N_alloc = N_alloc + allocsz;
               intpts(N_alloc,3) = 0;
               I(N_alloc) = 0;
               T(N_alloc) = 0;
               E(N_alloc) = 0;
            end
            idv = N_pts + (1:N);
            intpts(idv,:) = ipts(t,:);
            I(idv) = idv;
            T(idv) = fid(t);
            E(idv) = i;
            N_pts = N_pts + N;
         end
      end
   end
   T = T(1:N_pts);
   E = E(1:N_pts);
   I = I(1:N_pts);
   intpts = intpts(1:N_pts,:);
   tri2edge = sparse(T,E,I,size(fmdl.faces,1),size(rmdl.edges,1));
   tri2intpt = sparse(T,I,ones(size(I)),size(fmdl.faces,1),size(I,1));
   edge2intpt  = sparse(E,I,ones(size(I)),size(rmdl.edges,1),size(I,1));    
   
%-------------------------------------------------------------------------%
% Assign each rmdl node to the tet it is in (nodes on tet faces are counted
% mutltiple times)  
function rnode2tet = get_nodes_in_tets(fmdl,nodes, opt)
    
   [A,b] = tet_to_inequal(fmdl.nodes,fmdl.elems);
   progress_msg(.01);
   % This is split to decrease the memory footprint
   rnode2tet = (bsxfun(@minus, A(1:4:end,:)*nodes',b(1:4:end)) <= opt.tol_node2tet)';
   progress_msg(.21);
   for i = 2:4
      rnode2tet = rnode2tet & (bsxfun(@minus, A(i:4:end,:)*nodes',b(i:4:end)) <= opt.tol_node2tet)';
      progress_msg(.21 + (i-1)*.23);
   end

   % exclude coinciding nodes
   ex= bsxfun(@eq,nodes(:,1),fmdl.nodes(:,1)') & ...
       bsxfun(@eq,nodes(:,2),fmdl.nodes(:,2)') & ...
       bsxfun(@eq,nodes(:,3),fmdl.nodes(:,3)');
   progress_msg(.94);
   rnode2tet(any(ex,2),:) = 0;
   rnode2tet = sparse(rnode2tet);
   progress_msg(1);

%-------------------------------------------------------------------------%
% Prepare model
function fmdl = prepare_tet_mdl(fmdl)
   fmopt.elem2edge = true;
   fmopt.edge2elem = true;
   fmopt.face2elem = true;
   fmopt.node2elem = true;
   fmopt.normals   = true;
   fmopt.linear_reorder = false; % this is slow and not needed
   ll = eidors_msg('log_level',1);
   fmdl = fix_model(fmdl,fmopt);
   eidors_msg('log_level',ll);
   fmdl.node2elem = logical(fmdl.node2elem);
   nElem = size(fmdl.elems,1);
   nFace = size(fmdl.faces,1);
   fmdl.elem2face = sparse(repmat((1:nElem)',1,4),double(fmdl.elem2face),true,nElem,nFace);

%-------------------------------------------------------------------------%
% Center scale models
function[fmdl,rmdl] = center_scale_models(fmdl,rmdl, opt)
   ctr = mean([min(rmdl.nodes);max(rmdl.nodes)]);
   rmdl.nodes = bsxfun(@minus,rmdl.nodes,ctr);
   fmdl.nodes = bsxfun(@minus,fmdl.nodes,ctr);
   if isfield(opt,'do_not_scale') && opt.do_not_scale
      return
   end
   maxnode = min( max(abs(rmdl.nodes(:))), max(abs(fmdl.nodes(:))));
   scale = 1/maxnode;
   rmdl.nodes = scale*rmdl.nodes;
   fmdl.nodes = scale*fmdl.nodes;
   eidors_msg('@@ models scaled by %g', scale,2);

%-------------------------------------------------------------------------%
% Remove obviously non-overlapping parts of the models
function [fmdl,rmdl,fmdl_idx,rmdl_idx] = crop_models(fmdl,rmdl)
   f_min = min(fmdl.nodes);
   f_max = max(fmdl.nodes);
   r_min = min(rmdl.nodes);
   r_max = max(rmdl.nodes);
   
   % nodes outside the bounding box of the other model
   f_gt  = bsxfun(@gt, fmdl.nodes, r_max);
   f_lt  = bsxfun(@lt, fmdl.nodes, r_min);
   r_gt  = bsxfun(@gt, rmdl.nodes, f_max);
   r_lt  = bsxfun(@lt, rmdl.nodes, f_min);
   
   % elems outside the bounding box of the other model
   re_gt = any(reshape(all(reshape(r_gt(rmdl.elems',:),4,[])),[],3),2);
   re_lt = any(reshape(all(reshape(r_lt(rmdl.elems',:),4,[])),[],3),2);
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
   
   rmdl.elems = reshape(r_idx(r_n),[],4);
   fmdl.elems = reshape(f_idx(f_n),[],4);
   
   rmdl.nodes = rmdl.nodes(r_used_nodes,:);
   fmdl.nodes = fmdl.nodes(f_used_nodes,:);
   
   % for the benefit of any (debug) plots later on
   rmdl = rmfield(rmdl,'boundary');
   fmdl = rmfield(fmdl,'boundary');
    
%-------------------------------------------------------------------------%
% Parse option struct
 function opt = parse_opts(fmdl,rmdl, opt)

    
    if ~isfield(opt, 'tol_node2tet');
        opt.tol_node2tet = eps; % * max(rmdl_rng,fmdl_rng)^3;
    end
    if ~isfield(opt, 'tol_edge2edge')
        opt.tol_edge2edge = 2*sqrt(3)*eps(min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))));
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
   
   
%-------------------------------------------------------------------------%
% Perfom unit tests
function do_unit_test
   do_small_test;
   do_realistic_test;


function do_small_test
   fmdl = ng_mk_cyl_models([1 .5],[],[]);
   show_fem(fmdl)
   v = -.5:.1:.5;
   rmdl = mk_grid_model([],v,v,0:.1:1);
   hold on
   h = show_fem(rmdl);
   set(h,'edgecolor','b');
   hold off
   c2f = mk_tet_c2f(fmdl,rmdl);
   tc2f = c2f * rmdl.coarse2fine;
   vc2f = mk_grid_c2f(fmdl,rmdl);
   unit_test_cmp('mk_tet_c2f v mk_grid_c2f', tc2f,vc2f, 1e-15);


function do_realistic_test
   fmdl= ng_mk_cyl_models([2,2,.1],[16,1],[.1,0,.025]);
   xvec = [-1.5 -.5:.2:.5 1.5];
   yvec = [-1.6 -1:.2:1 1.6];
   zvec = 0:.25:2;
   rmdl = mk_grid_model([],xvec,yvec,zvec);
   tic
   opt.save_memory = 0;
   c2f_a = mk_grid_c2f(fmdl, rmdl,opt);
   t = toc;
   fprintf('Voxel: t=%f s\n',t);

   tic
   opt.save_memory = 0;
   c2f_b = mk_tet_c2f(fmdl, rmdl,opt);
   t = toc;
   fprintf('Tet: t=%f s\n',t);

   c2f_b = c2f_b * rmdl.coarse2fine;
   unit_test_cmp('mk_tet_c2f v mk_grid_c2f', c2f_b,c2f_a, 1e-5);

   tic
   c2f_n = mk_coarse_fine_mapping(fmdl,rmdl);
   t = toc;
   fprintf('Approximate: t=%f s\n',t);

