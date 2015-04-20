function [c2f] = mk_tet_c2f(fmdl, rmdl, opt)
%MK_TET_C2F

% (C) 2015 Bartlomiej Grychtol - all rights reserved by Swisstom AG
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<
if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'), do_unit_test; return; end
if nargin < 3
   opt = struct();
end

[fmdl,rmdl] = center_scale_models(fmdl,rmdl, opt);

opt = parse_opts(fmdl,rmdl, opt);


copt.fstr = 'mk_tet_c2f';

c2f = eidors_cache(@do_mk_tet_c2f,{fmdl,rmdl,opt},copt);

function c2f = do_mk_tet_c2f(fmdl,rmdl,opt)
   DEBUG = eidors_debug('query','mk_tet_c2f');
   
   c2f = sparse(0,0);
   fmdl = prepare_tet_mdl(fmdl);
   rmdl = prepare_tet_mdl(rmdl);
   progress_msg('Find c_edge2f_face intersections...')
   [intpts1, fface2redge, fface2intpt1, redge2intpt1] = ...
      edge2face_intersections(fmdl,rmdl);
   progress_msg(sprintf('Found %d', size(intpts1,1)), Inf);
   
   progress_msg('Find f_edge2c_face intersections...')
   [intpts2, rface2fedge, rface2intpt2, fedge2intpt2] = ...
      edge2face_intersections(rmdl,fmdl);
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
   rnode_in_ftet = (double(rmdl.node2elem') * rnode2ftet) == 4;
   progress_msg(sprintf('Found %d',nnz(rnode_in_ftet)), Inf);
   
   progress_msg('Find f_nodes in c_tets...');
   fnode2rtet = get_nodes_in_tets(rmdl,fmdl.nodes, opt);
   progress_msg(sprintf('Found %d', nnz(fnode2rtet)), Inf);

   progress_msg('Find f_elems in c_elems...')
   fnode_in_rtet = (double(fmdl.node2elem') * fnode2rtet) == 4;
   progress_msg(sprintf('Found %d',nnz(fnode_in_rtet)), Inf);
   
   progress_msg('Find total intersections...');
   e2e = double(rmdl.edge2elem');
   rtet2ftet =  double(rmdl.elem2face) * (rface2fedge>0) * fmdl.edge2elem ...
                 | e2e * (fface2redge>0)' * fmdl.elem2face' ...
                 | e2e * fedge2redge' * fmdl.edge2elem;
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
   
   for v = rtet_todo'
%         id = id+1;
%         if mod(id,mint)==0, progress_msg(id/lvox); end
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
            ok = false;
%             if size(pts,1) < 4 % test if edge lies on the plane of the vox
%                 % check for edges along the x y or z axis
%                 % this includes coplanar faces
%                 E = fmdl.edges(fmdl.elem2edge(tet_todo(t),:),:);
%                 P1 = fmdl.nodes(E(:,1),:);
%                 P2 = fmdl.nodes(E(:,2),:);
%                 % this test is sensitive, but not specific
%                 % it should also check if both pts come from the same edge and
%                 % that edge fullfils the condition
%                 D = P1-P2;
%                 ok = any(D(:) == 0); 
%             end 
            
            if ok, continue, end % otherwise convhulln will throw an error
            try
                % move points to origin (helps for small elements at
                % large coordinates
                ctr = mean(pts);
                pts = bsxfun(@minus,pts,ctr);
                scale = max(abs(pts(:)));
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
%                         % check for edges along the x y or z axis
%                         % this includes coplanar faces
%                         E = fmdl.edges(fmdl.elem2edge(tet_todo(t),:),:);
%                         P1 = fmdl.nodes(E(:,1),:);
%                         P2 = fmdl.nodes(E(:,2),:);
%                         % this test is sensitive, but not specific
%                         % it should also check if both pts come from the same edge and
%                         % that edge fullfils the condition
%                         D = P1-P2;
%                         ok = any(abs(D) < eps);
%                         % edge-edge intersections often appear to also
%                         % cross faces, there doesn't seem to be a good
%                         % specific way to catch that
%                         ok = ok | size(unique(pts,'rows'),1) < 4;
                end
                if ~ok
                    if DEBUG || eidors_debug('query','mk_tet_c2f:convhulln');
                        tet.nodes = fmdl.nodes;
                        vox.nodes = rmdl.nodes;
                        tet.type = 'fwd_model';
                        vox.type = 'fwd_model';
                        vox.elems = m.faces(logical(m.vox2face(v,:)),:);
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
    c2f = c2f + bsxfun(@times, sparse(rtet_in_ftet), get_elem_volume(rmdl))';
    
    % normalize to tet volume
    vol = get_elem_volume(fmdl);
    c2f = bsxfun(@rdivide,c2f,vol);

    % add tets contained in vox

    c2f = c2f + ftet_in_rtet;
%-------------------------------------------------------------------------%
% Calculate intersection points between faces and edges    
function [intpts, tri2edge, tri2intpt, edge2intpt] = edge2face_intersections(fmdl,rmdl,opt)
    
   intpts = [];
   T = []; E = []; I = [];
  
   P1 = rmdl.nodes(rmdl.edges(:,1),:);
   P12 = P1 - rmdl.nodes(rmdl.edges(:,2),:);
   N_edges = size(rmdl.edges,1);
   
   d = sum(fmdl.normals .* fmdl.nodes(fmdl.faces(:,1),:),2);
   
   num = -repmat(d,1,size(P1,1));
   for j = 1:3
      num = num + bsxfun(@times,fmdl.normals(:,j),P1(:,j)');
   end
   
   den = bsxfun(@times,fmdl.normals(:,j),P12(:,j)');
   for j = 2:3
      den = den + bsxfun(@times,fmdl.normals(:,j),P12(:,j)');
   end
   
   u = num ./ den;
   
   idx = u >= 0 & u <= 1;
   
   for i = find(any(idx))
      progress_msg(i/N_edges);
      % calculate the intersection points
      id = find(idx(:,i));
      ipts = bsxfun(@minus, P1(i,:), bsxfun(@times, u(id,i), P12(i,:)));
      t = point_in_triangle(ipts,fmdl.faces(id,:),fmdl.nodes,'match');
      if any(t)
         N = nnz(t);
         intpts = [intpts; ipts(t,:)];
         I = [I; (1:N)' + size(I,1)];
         T = [T; id(t)];
         E = [E; i*ones(N,1)];
      end
   end
   
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
   progress_msg(1);

%-------------------------------------------------------------------------%
% Prepare matrices for the voxel model
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


function do_small_test
   fmdl = ng_mk_cyl_models([1 .5],[],[]);
   show_fem(fmdl)
   v = -.5:.1:.5;
   rmdl = mk_grid_model([],v,v,0:.1:1);
   hold on
   h = show_fem(rmdl);
   set(h,'edgecolor','b');
   mk_tet_c2f(fmdl,rmdl);
