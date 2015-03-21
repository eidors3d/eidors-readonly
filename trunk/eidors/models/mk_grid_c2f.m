function [c2f, m] = mk_grid_c2f(fmdl, rmdl, opt)
%MK_GRID_C2F - calculate a coarse2fine mapping for grid coarse models.
% C2F = MK_GRID_C2F(FMDL,RMDL) returns in C2F the fraction of volume of
% each element of the fine (tet-based) model contained in each element of
% the coarse (vox-based) model.
% Uses CONVHULLN to calculate the volume defined by a set of intersection
% points between individual tet and vox elements.
%
% C2F = MK_GRID_C2F(FMDL,RMDL,OPT) allows specifying options.
% 
% Inputs:
%   FMDL - an EIDORS (tet-based) forward model
%   RMDL - a grid model, as returned by MK_GRID_MODEL
%   OPT  - an option structure with the following fields and defaults:
%      .do_not_scale  - set to true to prevent scaling the models to unit
%                       cube before any calculations, including thresholds.
%                       Default: false
%      .tol_node2tet  - tolerance for determinant <= 0 in testing for
%                       points inside tets. Default: eps
%      .tol_edge2edge - maximum distance between "intersecting" edges
%                       Default: 2*sqrt(3)*eps(a), where a is
%                       min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))
%      .tol_edge2tri  - minimum value of a barycentric coordinate to 
%                       decide a point is lying inside a triangle and not
%                       on its edge. Default: eps
%      .save_memory   - modifies function behavior to sacrifice time to
%                       decrease memory footprint, useful for large
%                       problems. Must be an integer between 0 and 3
%                          0  -  calculate all at once (default)
%                          1  -  calculate one xy voxel plane at a time
%                          2  -  calculate one y voxel row at a time
%                          3  -  calculate each voxel separately
%
% [C2F M] = MK_GRID_C2F(...) also returns a struct with useful fields
% charactrising the vox model.
%
% Set eidors_msg 'log level' < 2 to supress output to command line.
%
% Examples:
%     fmdl = ng_mk_cyl_models([2,2,.2],[],[]);
%     rmdl = mk_grid_model([],-2:2,-2:2,0:2);
%     c2f  = mk_grid_c2f(fmdl,rmdl);
%     h = show_fem(fmdl); set(h,'LineWidth',0.1)
%     hold on
%     h = show_fem(rmdl); set(h,'EdgeColor','b','LineWidth',2);
%     hold off
%
% See also MK_GRID_MODEL, FIND_EDGE2EDGE_INTERSECTIONS, CONVHULLN
%          MK_COARSE_FINE_MAPPING, EIDORS_MSG

% (C) 2015 Bartlomiej Grychtol - all rights reserved by Swisstom AG
% License: GPL version 2 or 3
% $Id$

% >> SWISSTOM CONTRIBUTION <<


% if input is 'UNIT_TEST', run tests
if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'), do_unit_test; return; end
if ischar(fmdl) && strcmp(fmdl,'PROFILE'), do_small_test; return; end
if nargin < 3
    opt = struct();
end

[fmdl,rmdl] = center_scale_models(fmdl,rmdl, opt);

opt = parse_opts(fmdl,rmdl, opt);

copt.cache_obj = {fmdl.nodes,
                  fmdl.elems,
                  rmdl.nodes,
                  rmfield(opt,'save_memory')};
               
copt.fstr = 'mk_grid_c2f';

[c2f, m] = eidors_cache(@do_mk_grid_c2f,{fmdl,rmdl,opt},copt);




function [c2f, m]= do_mk_grid_c2f(fmdl0,rmdl0,opt0)
    eidors_msg('@@@',2);
    
    m = [];
    c2f = sparse(opt0.nTet,opt0.nVox);
    
    logmsg(' Prepare tet model... ');
    fmdl0 = prepare_fmdl(fmdl0);
    logmsg('Done\n');
    
    if opt0.save_memory == 0
       relem_idx = ones(opt0.nVox,1);
       [fmdl, rmdl, opt, felem_idx] = crop_models(fmdl0,rmdl0,opt0,relem_idx);
       if any(felem_idx)
          [tmp, m] = separable_calculations(fmdl,rmdl0,opt);
          c2f = combine_c2f(c2f, tmp,felem_idx,relem_idx,opt,opt0);
       end
    else % save_memory > 0
       logmsg(' Saving memory mode level %d\n',opt0.save_memory);
       logmsg(' Progress: ');
       progmsg(0,false);
       progress = 0;
       opt = opt0;
       m = prepare_vox_mdl(rmdl0,opt0);
       for z = 1:opt0.Zsz
          opt.Zsz = 1;
          opt.zvec = opt0.zvec(z:z+1);
          opt.xplane = opt0.xplane(z,:);
          opt.yplane = opt0.yplane(z,:);
          relem_idx = false(opt0.Xsz*opt0.Ysz*opt0.Zsz,1);
          relem_idx((z-1)*opt0.Xsz*opt0.Ysz + (1:opt0.Xsz*opt0.Ysz)) = true;
          if opt0.save_memory > 1
             for y = 1:opt0.Ysz
                opt.Ysz = 1;
                opt.yvec = opt0.yvec(y:y+1);
                opt.zplane = opt0.zplane(y,:);
                opt.xplane = opt0.xplane(z,y);
                opt.zstep = opt.ystep*(opt.Ysz+1);
                relem_idx_y = false(opt0.Xsz*opt0.Ysz,1);
                relem_idx_y((y-1)*opt0.Xsz + (1:opt0.Xsz)) = true;
                relem_idx_y = relem_idx & repmat(relem_idx_y,opt0.Zsz,1);
                if opt0.save_memory > 2
                   for x = 1:opt0.Xsz
                      opt.Xsz = 1;
                      opt.xvec = opt0.xvec(x:x+1);
                      opt.yplane = opt0.yplane(z,x);
                      opt.zplane = opt0.zplane(y,x);
                      opt.ystep = opt.xstep*opt.Xsz+1;
                      opt.zstep = opt.ystep*(opt.Ysz+1);
                      relem_idx_x = false(opt0.Xsz,1);
                      relem_idx_x(x) = true;
                      relem_idx_x = relem_idx_y & repmat(relem_idx_x,opt0.Ysz*opt0.Zsz,1);
                      [fmdl, rmdl, opt, felem_idx] = crop_models(fmdl0,rmdl0,opt,relem_idx_x);
                      if any(felem_idx)
                         eidors_msg('log_level',eidors_msg('log_level')-1);
                         tmp = separable_calculations(fmdl,rmdl,opt);
                         eidors_msg('log_level',eidors_msg('log_level')+1);
                         c2f = combine_c2f(c2f, tmp,felem_idx,relem_idx_x,opt,opt0);
                      end
                      progress = progress + 1 / opt0.nVox;
                      progmsg(progress);
                   end
                else
                   [fmdl, rmdl, opt, felem_idx] = crop_models(fmdl0,rmdl0,opt,relem_idx_y);
                   if any(felem_idx)
                      eidors_msg('log_level',eidors_msg('log_level')-1);
                      tmp = separable_calculations(fmdl,rmdl,opt);
                      eidors_msg('log_level',eidors_msg('log_level')+1);
                      c2f = combine_c2f(c2f, tmp,felem_idx,relem_idx_y,opt,opt0);
                   end
                   progress = progress + 1 / (opt0.Xsz*opt0.Ysz);
                   progmsg(progress);
                end
             end

          else
             [fmdl, rmdl, opt, felem_idx] = crop_models(fmdl0,rmdl0,opt,relem_idx);
             if any(felem_idx)
                eidors_msg('log_level',eidors_msg('log_level')-1);
                tmp = separable_calculations(fmdl,rmdl,opt);
                eidors_msg('log_level',eidors_msg('log_level')+1);
                c2f = combine_c2f(c2f, tmp,felem_idx,relem_idx,opt,opt0);
             end
             progress = progress + 1 / opt0.Zsz;
             progmsg(progress);
          end
          
          
       end
       logmsg('\b\b\b\bDone\n');
    end
    
function c2f = combine_c2f(c2f, tmp,felem_idx,relem_idx,opt,opt0)
    F = repmat(find(felem_idx), 1, opt.nVox);
    R = repmat(find(relem_idx)', opt.nTet, 1);
    c2f = c2f + sparse(F, R, tmp, opt0.nTet, opt0.nVox);
    
function [c2f, m] = separable_calculations(fmdl,rmdl,opt)
    DEBUG = eidors_debug('query','mk_grid_c2f');
   
    logmsg(' Prepare vox model... ');
    m = prepare_vox_mdl(rmdl,opt);
    logmsg('Done\n');
    
    try
    
    % tet edge v. vox face
    logmsg(' Find tet_edge2vox_face intersections... ')
    [intpts1, rec2tedge, rec2intpt1, tedge2intpt1] = get_voxel_intersection_points(fmdl,m.faces,opt);
    logmsg('\b\b\b\bFound %d\n', size(intpts1,1));
    
    % vox edge v. tet face
    logmsg(' Find vox_edge2tet_face intersections... ')
    [intpts2, tri2vedge, tri2intpt2, vedge2intpt2] = get_tet_intersection_points(fmdl,m,opt);
    logmsg('\b\b\b\bFound %d\n', size(intpts2,1));
    
    % Note: Rather than calculating edge2edge intersections, one could
    % include them in one or both of the previous tests. However, it is
    % then difficult to guarantee that the numerically border-line cases
    % get assigned to all the vox and tets concerned
    
    % vox edge v. tet edge
    logmsg(' Find tet_edge2vox_edge intersections... ')
    [intpts3, tedge2vedge, tedge2intpt3, vedge2intpt3] =find_edge2edge_intersections(fmdl.edges,fmdl.nodes,m.edges,rmdl.nodes, opt.tol_edge2edge);
    logmsg('\b\b\b\b\bFound %d\n', size(intpts3,1));
    
    % tet node in vox
    logmsg(' Find tet_nodes in voxels... ')
    [fnode2vox] = get_nodes_in_voxels(fmdl,rmdl);
    logmsg('Found %d\n', nnz(fnode2vox));
    % vox node in tet
    logmsg(' Find vox_nodes in tets... ');
    vnode2tet = get_nodes_in_tets(fmdl,rmdl, opt);
    logmsg('\b\b\b\bFound %d\n', nnz(vnode2tet));
    
    % vox contained in tet
    vox_in_tet = (m.node2vox' * vnode2tet) == 8;
    logmsg(' Find vox contained in tet... Found %d\n',nnz(vox_in_tet));
    
    % tet contained in vox
    tet_in_vox = (double(fmdl.node2elem') * fnode2vox) == 4;
    logmsg(' Find tets contained in vox... Found %d\n',nnz(tet_in_vox));
    
    
    % tets and vox that intersect (excludes complete inclusion)
    logmsg(' Find total vox v. tet intersections... ');
    vox2intTet =   m.vox2face * (rec2tedge>0) * fmdl.edge2elem ...
                 | m.vox2edge * (tri2vedge>0)' * fmdl.elem2face';% ...
%                  | m.vox2edge * (tedge2vedge>0)'* fmdl.edge2elem;
    logmsg('Found %d\n',nnz(vox2intTet));
    
    catch err
       if (strcmp(err.identifier,'MATLAB:nomem'))
          msg = sprintf('%s', ...
             'Matlab ran out of memory. Consider setting opt.save_memory > 0.\n', ...
             'See HELP MK_GRID_C2F for details.');
          error('EIDORS:mk_grid_c2f:nomem',msg);
       else
          rethrow(err);
       end
    end
    
    logmsg(' Calculate intersection volumes... ');
    % sparse logical multiplication doesn't exist
    vox2intpt1 = logical(m.vox2face*rec2intpt1)'; 
    tet2intpt1 = logical(fmdl.edge2elem'*tedge2intpt1)';

    tet2intpt2 = logical(fmdl.elem2face*tri2intpt2)';
    vox2intpt2 = logical(m.vox2edge*vedge2intpt2)';

    tet2intpt3 = logical(fmdl.edge2elem'*tedge2intpt3)';
    vox2intpt3 = logical(m.vox2edge*vedge2intpt3)';
    
    vox_todo = find(sum(vox2intTet,2)>0);
    C = []; F = []; V = [];
    progmsg(0,false);
    id = 0; lvox = length(vox_todo);
    mint = ceil(lvox/100);
    for v = vox_todo'
        id = id+1;
        if mod(id,mint)==0, progmsg(id/lvox); end
        tet_todo = find(vox2intTet(v,:));
        common_intpts1 = bsxfun(@and,vox2intpt1(:,v), tet2intpt1(:,tet_todo));
        common_intpts2 = bsxfun(@and,vox2intpt2(:,v), tet2intpt2(:,tet_todo));
        common_intpts3 = bsxfun(@and,vox2intpt3(:,v), tet2intpt3(:,tet_todo));
        tet_nodes     = bsxfun(@and,fnode2vox(:,v), fmdl.node2elem(:,tet_todo));
        vox_nodes     = bsxfun(@and,vnode2tet(:,tet_todo), m.node2vox(:,v));
        C = [C; v*ones(numel(tet_todo),1)];
        F = [F; tet_todo'];
        last_v = numel(V);
        V = [V; zeros(numel(tet_todo),1)]; % pre-allocate

        for t = 1:numel(tet_todo)
            pts = [ intpts1(common_intpts1(:,t),:);
                    intpts2(common_intpts2(:,t),:);
                    intpts3(common_intpts3(:,t),:);
                    fmdl.nodes(tet_nodes(:,t),:);
                    rmdl.nodes(vox_nodes(:,t),:)];
            last_v = last_v + 1;
            ok = false;
            if size(pts,1) < 4 % test if edge lies on the plane of the vox
                % check for edges along the x y or z axis
                % this includes coplanar faces
                E = fmdl.edges(fmdl.elem2edge(tet_todo(t),:),:);
                P1 = fmdl.nodes(E(:,1),:);
                P2 = fmdl.nodes(E(:,2),:);
                % this test is sensitive, but not specific
                % it should also check if both pts come from the same edge and
                % that edge fullfils the condition
                D = P1-P2;
                ok = any(D(:) == 0); 
            end 
            
            if ok, continue, end % otherwise convhulln will throw an error
            try
                % move points to origin (helps for small elements at
                % large coordinates
                ctr = mean(pts);
                pts = bsxfun(@minus,pts,ctr);
                scale = max(abs(pts(:)));
                % scale largest coordinate to 1 (helps with precision)
                pts = pts ./ scale;
                % force thourugh search for initinal simplex and 
                % supress precision warnings
                [K, V(last_v)] = convhulln(pts,{'Qt Pp Qs'});
                V(last_v) = V(last_v) * scale^3; % undo scaling
            catch err
                ok = false;
                switch err.identifier
                    case {'MATLAB:qhullmx:DegenerateData', 'MATLAB:qhullmx:UndefinedError'}
                        % check for edges along the x y or z axis
                        % this includes coplanar faces
                        E = fmdl.edges(fmdl.elem2edge(tet_todo(t),:),:);
                        P1 = fmdl.nodes(E(:,1),:);
                        P2 = fmdl.nodes(E(:,2),:);
                        % this test is sensitive, but not specific
                        % it should also check if both pts come from the same edge and
                        % that edge fullfils the condition
                        D = P1-P2;
                        ok = any(abs(D) < eps);
                        % edge-edge intersections often appear to also
                        % cross faces, there doesn't seem to be a good
                        % specific way to catch that
                        ok = ok | size(unique(pts,'rows'),1) < 4;
                end
                if ~ok
                    if DEBUG || eidors_debug('query','mk_grid_c2f:convhulln');
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
                            'Enable eidors_debug on mk_grid_c2f and re-run to see a debug plot'],0);
                        rethrow(err);
                    end
                end
            end
        end
    end
    logmsg('\b\b\b\bDone\n');
    c2f = sparse(F,C,V,opt.nTet,opt.nVox);
    
    % add vox contained in tet
    c2f = c2f + bsxfun(@times, sparse(vox_in_tet), m.volume)';
    
    % normalize to tet volume
    vol = get_elem_volume(fmdl);
    c2f = bsxfun(@rdivide,c2f,vol);

    % add tets contained in vox

    c2f = c2f + tet_in_vox;
    
function [fmdl, rmdl, opt, felem_idx] = crop_models(fmdl0,rmdl0,opt, relem_idx)
   
   fmdl = fmdl0;
   rmdl = rmdl0;
   
   opt.nVox = opt.Xsz*opt.Ysz*opt.Zsz;
   idx = rmdl0.coarse2fine * relem_idx > 0;
   rmdl.elems = rmdl0.elems(idx,:);
   [node_idx, m,n] = unique(rmdl.elems(:));
   rmdl.elems = reshape(n,size(rmdl.elems));
   rmdl.nodes = rmdl0.nodes(node_idx,:);
   
   fnode_above = fmdl0.nodes(:,3) > opt.zvec(end);
   felem_above = all(fnode_above(fmdl0.elems),2);
   
   fnode_below = fmdl0.nodes(:,3) < opt.zvec(1);
   felem_below = all(fnode_below(fmdl0.elems),2);
   
   felem_idx   = ~(felem_below | felem_above);
   fmdl.elems = fmdl0.elems(felem_idx,:);
   fnode_idx = unique(fmdl.elems(:));
   fnode_idx_map = zeros(size(fmdl0.nodes,1),1);
   fnode_idx_map(fnode_idx) = 1:length(fnode_idx);
   fmdl.elems = fnode_idx_map(fmdl.elems);
   fmdl.edges = fnode_idx_map(fmdl0.edges);
   fmdl.nodes = fmdl0.nodes(fnode_idx,:);
   fmdl.node2elem = fmdl0.node2elem(fnode_idx,felem_idx);
   
   
   fface_idx = sum(fmdl0.elem2face(felem_idx,:))>0;
   fmdl.elem2face = fmdl0.elem2face(felem_idx,fface_idx);
   felem_idx_map = felem_idx;
   felem_idx_map(felem_idx) = 1:nnz(felem_idx);
   felem_idx_map = [0; felem_idx_map];
   fmdl.face2elem = fmdl0.face2elem(fface_idx,:);
   fmdl.face2elem = felem_idx_map(fmdl.face2elem + 1);
   fmdl.normals = fmdl0.normals(fface_idx,:);
   fmdl.faces = fmdl0.faces(fface_idx,:);
   fmdl.faces = fnode_idx_map(fmdl.faces);
   
   fedge_idx = unique(fmdl0.elem2edge(felem_idx,:));
   fedge_idx_map = zeros(size(fmdl0.edges,1),1);
   fedge_idx_map(fedge_idx) = 1:length(fedge_idx);
   fmdl.elem2edge = fedge_idx_map(fmdl0.elem2edge(felem_idx,:));
   fmdl.edge2elem = fmdl0.edge2elem(fedge_idx,felem_idx);
   fmdl.edges = fmdl.edges(fedge_idx,:);
   
   opt.nTet = size(fmdl.elems,1);

%-------------------------------------------------------------------------%
% Prepare matrices for the voxel model
function fmdl = prepare_fmdl(fmdl)
    fmopt.elem2edge = true;
    fmopt.edge2elem = true;
    fmopt.face2elem = true;
    fmopt.node2elem = true;
    fmopt.normals   = true;
    fmopt.linear_reorder = false; % this is slow and not needed
    fmdl = fix_model(fmdl,fmopt);
    fmdl.node2elem = logical(fmdl.node2elem);
    nElem = size(fmdl.elems,1);
    nFace = size(fmdl.faces,1);
    fmdl.elem2face = sparse(repmat((1:nElem)',1,4),double(fmdl.elem2face),true,nElem,nFace);


%-------------------------------------------------------------------------%
% Prepare matrices for the voxel model
function m = prepare_vox_mdl(rmdl,opt)

    DEBUG = eidors_debug('query','mk_grid_c2f');

    [voxels, m.node2vox] = mk_voxels(opt);
    if DEBUG || eidors_debug('query','mk_grid_c2f:mk_voxels')
        show_voxels(rmdl,voxels); title('mk\_voxels');
    end

    m.faces = mk_faces(voxels,opt);
    if DEBUG || eidors_debug('query','mk_grid_c2f:mk_faces')
        test_faces(rmdl,voxels(m.faces,:),opt); title('mk\_faces');
    end

    m.vox2face = mk_vox2face(opt);
    if DEBUG || eidors_debug('query','mk_grid_c2f:mk_vox2face')
        % the numbers shown are useless, just check if all faces are present
        show_voxels(rmdl,m.faces(any(m.vox2face),:));title('mk\_vox2face');
    end
    m.edges = mk_edges(voxels,opt);

    m.edge_length = rmdl.nodes(m.edges(:,1),:) - rmdl.nodes(m.edges(:,2),:);
    m.edge_length = sqrt(sum(m.edge_length.^2,2));
    
    [m.vox2edge, m.volume]= mk_vox2edge(m,opt);


%-------------------------------------------------------------------------%
% Assign each rmdl node to the tet it is in (nodes on tet faces are counted
% mutltiple times)  
function rnode2tet = get_nodes_in_tets(fmdl,rmdl, opt)
    
    [A,b] = tet_to_inequal(fmdl.nodes,fmdl.elems);
    progmsg(.01,false);
    % This is split to decrease the memory footprint
    rnode2tet = (bsxfun(@minus, A(1:4:end,:)*rmdl.nodes',b(1:4:end)) <= opt.tol_node2tet)';
    progmsg(.21);
    for i = 2:4 
        rnode2tet = rnode2tet & (bsxfun(@minus, A(i:4:end,:)*rmdl.nodes',b(i:4:end)) <= opt.tol_node2tet)'; 
        progmsg(.21 + (i-1)*.23);
    end
    
    % exclude coinciding nodes
    ex= bsxfun(@eq,rmdl.nodes(:,1),fmdl.nodes(:,1)') & ...
        bsxfun(@eq,rmdl.nodes(:,2),fmdl.nodes(:,2)') & ...
        bsxfun(@eq,rmdl.nodes(:,3),fmdl.nodes(:,3)');
    progmsg(.94);
    rnode2tet(any(ex,2),:) = 0;
    progmsg(1);
    

%-------------------------------------------------------------------------%
% Assign each fmdl node to the vox it is in (nodes on vox faces are counted
% mutltiple times)  
function [insnode] = get_nodes_in_voxels(fmdl,rmdl)

    E = reshape(rmdl.elems',4*6,[])';
    E = E(:,[1 2 3 4 8 12 16 23]);

    NE = size(E,1);
    xnodes = reshape(rmdl.nodes(E,1),NE,[]);
    ynodes = reshape(rmdl.nodes(E,2),NE,[]);
    znodes = reshape(rmdl.nodes(E,3),NE,[]);
    minx = min(xnodes,[],2);
    maxx = max(xnodes,[],2);
    miny = min(ynodes,[],2);
    maxy = max(ynodes,[],2);
    minz = min(znodes,[],2);
    maxz = max(znodes,[],2);

    leftof  = bsxfun(@lt, fmdl.nodes(:,1), minx');
    rightof = bsxfun(@gt, fmdl.nodes(:,1), maxx');
    infront = bsxfun(@lt, fmdl.nodes(:,2), miny');
    behind  = bsxfun(@gt, fmdl.nodes(:,2), maxy');
    below   = bsxfun(@lt, fmdl.nodes(:,3), minz');
    above   = bsxfun(@gt, fmdl.nodes(:,3), maxz');

    outnode = leftof | rightof | behind | infront | below | above;
    insnode = sparse(~outnode);

%-------------------------------------------------------------------------%
% Calculate intersection points between vox faces and tet edges    
function [intpts, face2edge, face2intpt, edge2intpt] = get_voxel_intersection_points(fmdl,faces,opt)
edges = fmdl.edges;
nodes = fmdl.nodes;
dir   = nodes(edges(:,2),:) - nodes(edges(:,1),:);
intpts = [];
F = []; E = []; I = [];

SZ = [opt.Xsz, opt.Ysz, opt.Zsz];
VEC = {opt.xvec, opt.yvec, opt.zvec};
% show_voxels(rmdl, voxels(faces,:))

%  mdl = rmdl;
%  mdl.elems = voxels(faces,:);
%  img = mk_image(mdl,1);
todo = sum(SZ)+3;
id = 0;
step = ceil(todo/100);
progmsg(0,false);
for d = 1:3
    for x = 0:SZ(d)
        id = id+1;
        if mod(id,step)==0, progmsg(id/todo); end
        % plane-edge intersection
        t = (VEC{d}(x+1) - nodes(edges(:,1),d)) ./ dir(:,d);
        crossing = t>0 & t<1;
        crossed  = find(crossing);
        if isempty(crossed), continue, end;
        % intersection points
        tmp   = nodes(edges(crossing,1),:) + bsxfun(@times,t(crossing),dir(crossing,:));
        face_coord = zeros(length(crossed),3);
        rmv = false(length(crossed),1);
        for i = 1:3
            if i == d
               face_coord(:,i) = x+1;
            else
                face_coord(:,i) = sum(bsxfun(@times,diff(bsxfun(@lt, tmp(:,i), VEC{i}'),1,2),(1:SZ(i))),2);
                rmv = rmv | face_coord(:,i) == 0;
            end
        end
        % also reject intersections with voxel nodes and edges
        same_x = any(bsxfun(@eq,tmp(:,1),opt.xvec'),2);
        same_y = any(bsxfun(@eq,tmp(:,2),opt.yvec'),2);
        same_z = any(bsxfun(@eq,tmp(:,3),opt.zvec'),2);
        rmv = rmv | (same_x + same_y + same_z) > 1;
        if nnz(rmv) == numel(rmv), continue, end
        tmp(rmv,:) = []; 
        face_coord(rmv,:) = [];
        crossed(rmv) = [];
        face_coord = face_coord - 1;
        face_idx = d + face_coord(:,1)*opt.xstep + face_coord(:,2)*opt.ystep + face_coord(:,3)*opt.zstep;
        I = [I; (1:size(tmp,1))' + size(I,1)];
        F = [F; face_idx];
        E = [E; crossed];
        intpts = [intpts; tmp];
%         unique(face_idx)
%         img.elem_data(:) = 0;
%         img.elem_data(unique(face_idx)) = 1;
%         show_fem(img,[0 0 1]);
        %     hold on
        %     x1 =fmdl.nodes(fmdl.edges(crossing,1),1);
        %     x2 =fmdl.nodes(fmdl.edges(crossing,2),1);
        %     y1 =fmdl.nodes(fmdl.edges(crossing,1),2);
        %     y2 =fmdl.nodes(fmdl.edges(crossing,2),2);
        %     z1 =fmdl.nodes(fmdl.edges(crossing,1),3);
        %     z2 =fmdl.nodes(fmdl.edges(crossing,2),3);
        %     plot3([x1 x2]',[y1 y2]',[z1 z2]','b');
        %     plot3(intpts(:,1),intpts(:,2),intpts(:,3),'ro','MarkerFaceColor','r');
        %     hold off
        %     axis auto
        %     view(2)
        %     pause
    end
end
face2edge = sparse(F,E,I,size(faces,1),size(edges,1));
face2intpt = sparse(F,I,ones(size(I)),size(faces,1),size(I,1));
edge2intpt  = sparse(E,I,ones(size(I)),size(edges,1),size(I,1));

%-------------------------------------------------------------------------%
% Calculate intersection points between tet faces and vox edges    
function [intpts, tri2edge, tri2intpt, edge2intpt] = get_tet_intersection_points(fmdl,m,opt)
    
    intpts = [];
    T = []; E = []; I = [];
    Xsz = opt.Xsz; Ysz = opt.Ysz; Zsz = opt.Zsz;
    SZ = [opt.Xsz, opt.Ysz, opt.Zsz];
    VEC = {opt.xvec, opt.yvec, opt.zvec};
    STEP(1) = opt.xstep; 
    STEP(2) = STEP(1)*(Xsz+1);
    STEP(3) = STEP(2)*(Ysz+1);
    
    d = sum(fmdl.normals .* fmdl.nodes(fmdl.faces(:,1),:),2);
    progmsg(0,false);
    line_axis = [3 1 2];
    for i = 1:3
        progmsg((i-1)/3);
        % define edge lines
        idx = 1:3;
        op = line_axis(i);
        idx(op) = [];
        
        pts = [repmat(VEC{idx(1)},SZ(idx(2))+1,1) kron(VEC{idx(2)},ones(SZ(idx(1))+1,1)) ];
        pt_idx = uint32(repmat((0:SZ(idx(1)))',SZ(idx(2))+1,1)*STEP(idx(1)) ...
                + kron((0:SZ(idx(2)))', ones(SZ(idx(1))+1,1))*STEP(idx(2)));
        
        % project on faces
        % plane equation is ax+by+cz+d = 0, where d = -(ax0 + by0 + cz0)
        z = repmat(d,1,size(pts,1));
        for j = 1:2
            z = z - bsxfun(@times,fmdl.normals(:,idx(j)),pts(:,j)');
        end
        z = z ./ repmat(fmdl.normals(:,op),1,size(pts,1));
        in = point_in_triangle(pts,fmdl.faces,fmdl.nodes(:,idx),-opt.tol_edge2tri)';

        
        % reject voxel nodes
        v = VEC{op}';
        in = in & ~reshape(any(bsxfun(@eq,z(:),v),2),size(in));
        M = bsxfun(@lt, z(:),v);
        M = xor(M(:,1:end-1), M(:,2:end));
        edge_num = reshape(uint32(sum(bsxfun(@times,uint32(M),uint32(1:SZ(op))),2)), size(z));
        in = in & edge_num;
        if nnz(in) == 0, continue, end
        edge_num(~in) = 0;

        edge_num(in) = (edge_num(in)-1) * STEP(op);
        edge_idx = edge_num + bsxfun(@times,uint32(in), uint32(i) + pt_idx');

        [t, p] = find(in);
        tmp = zeros(length(p),3);
        tmp(:,idx) = pts(p,1:2);
        tmp(:,op) = z(in);
        
        I = [I; (1:size(tmp,1))' + size(I,1)];
        T = [T; t];
        E = [E; edge_idx(in)];
        intpts = [intpts; tmp];        
    end
    
    tri2edge = sparse(T,E,I,size(fmdl.faces,1),size(m.edges,1));
    tri2intpt = sparse(T,I,ones(size(I)),size(fmdl.faces,1),size(I,1));
    edge2intpt  = sparse(E,I,ones(size(I)),size(m.edges,1),size(I,1));

   
%-------------------------------------------------------------------------%
% Make voxels
function [voxels, node2vox] = mk_voxels(opt)
    Xsz = opt.Xsz; Ysz = opt.Ysz; Zsz = opt.Zsz;
    Xp = Xsz+1; Yp = Ysz+1; Zp = Zsz+1; % number of planes

    up = Xp*Yp;
    vox = [ 1       1+up    1+up+Xp     1+Xp;
            1       2       2+up        1+up;
            1       1+Xp    2+Xp        2;
            2       2+up    2+up+Xp     2+Xp;
            1+Xp    2+Xp    2+up+Xp     1+up+Xp;
            1+up    1+Xp+up 2+Xp+up     2+up];

    voxrow   = bsxfun(@plus, repmat(vox,     Xsz,1), kron(0:Xsz-1          ,ones(1,6))');
    voxplane = bsxfun(@plus, repmat(voxrow,  Ysz,1), kron(Xp*(0:Ysz-1) , ones(1, 6*Xsz))');
    voxels   = bsxfun(@plus, repmat(voxplane,Zsz,1), kron(up*(0:Zsz-1) , ones(1, 6*Xsz*Ysz))');
    nVox     = Xsz * Ysz * Zsz;
    node2vox = sparse(voxels(1:3:end,:),kron((1:nVox)',ones(2,4)),1);
    
%-------------------------------------------------------------------------%
% Make faces
function faces = mk_faces(voxels,opt)
    Xsz = opt.Xsz; Ysz = opt.Ysz; Zsz = opt.Zsz;
    
    facerow     = [nonzeros(bsxfun(@plus,(1:3)',0:6:(Xsz-1)*6)); (Xsz-1)*6+4];
    faceplane   = bsxfun(@plus,facerow,0:6*Xsz:6*Xsz*(Ysz-1));
    % cap duplicates faces in order to preserve nice indexing
    cap         = kron(ones(3,1),faceplane(2:3:end,end)')+3;
    faceplane   = [faceplane,[ cap(:); cap(end)]];
    faces       = bsxfun(@plus, faceplane(:), 0:6*Xsz*Ysz:6*Xsz*Ysz*(Zsz-1));
    % cap duplicates faces in order to preserve nice indexing
    cap         = reshape(kron(ones(3,1), 3:3:3*Xsz*Ysz'),3*Xsz,[]);
    cap         = bsxfun(@plus, cap, 0:Ysz-1);
    cap(end+1,:)= cap(end,:);
    faces       = [faces(:); faces(cap(:),end)+3];
    faces       = voxels(faces,:);
    
%-------------------------------------------------------------------------%
% Make edges
function edges = mk_edges(voxels, opt)
    Xsz = opt.Xsz; Ysz = opt.Ysz; Zsz = opt.Zsz;
    
    edgerow     = voxels([nonzeros(bsxfun(@plus,(1:3)',0:6:(Xsz-1)*6)); (Xsz-1)*6+4],1:2);
    edgerow(end+1,:) = edgerow(end,:);
    edgerow(end+1,:) = voxels((Xsz-1)*6+4, [1 4]);
    edgeplane = repmat(edgerow,Ysz+1,1) + kron((Xsz+1)*(0:Ysz)', ones(size(edgerow)));
    % replace ficticious edges with repetitions;
    edgeplane((size(edgerow,1)*Ysz +3):3:end,:) = edgeplane((size(edgerow,1)*Ysz +2):3:end,:);

    up = (Xsz+1)*(Ysz+1);
    edges = repmat(edgeplane,Zsz+1,1) + kron((0:Zsz)'*up, ones(size(edgeplane)));
    
    % replace ficticious edges with repetitions;
    start = (size(edgeplane,1)*Zsz);
    edges( (start +1):3:end,:) = edges( (start +2):3:end,:);
    start = (size(edgeplane,1)*Zsz) + 3*Xsz;
    step  = size(edgerow,1);
    edges( (start +1):step:end,:) = edges( (start +3):step:end,:);
    edges( (start +2):step:end,:) = edges( (start +3):step:end,:);
    edges( end-2:end,:) = edges(end-5:end-3,:);

 %-------------------------------------------------------------------------%
% Make mapping between voxels and faces  
function vox2face = mk_vox2face(opt)
    Xsz = opt.Xsz; Ysz = opt.Ysz; Zsz = opt.Zsz;
    xstep = opt.xstep; ystep = opt.ystep; zstep = opt.zstep;
    
    vox = [1 2 3 1+xstep 2+ystep 3+zstep];
    voxrow = bsxfun(@plus,(0:Xsz-1)'*xstep,vox)';
    voxplane = bsxfun(@plus,(0:Ysz-1)*ystep,voxrow(:));
    voxvol   = bsxfun(@plus,(0:Zsz-1)*zstep,voxplane(:));
    voxvol   = reshape(voxvol, 6, [])';
    I = kron((1:size(voxvol,1))',ones(1,6));
    nFaces  = Zsz*(Ysz+1)*(3*Xsz+1) + Ysz*(3*Xsz+1);
    nVox    = Xsz*Ysz*Zsz;
    vox2face = sparse(I,voxvol,ones(size(I)),nVox,nFaces);   
    
%-------------------------------------------------------------------------%
% Make mapping between voxels and edges  
function [vox2edge, vol] = mk_vox2edge(m,opt)
    Xsz = opt.Xsz; Ysz = opt.Ysz; Zsz = opt.Zsz;
    xstep = opt.xstep; ystep = xstep*(Xsz+1); zstep = ystep*(Ysz+1);
    
    vox = [1 2 3 1+xstep 3+xstep 1+ystep 2+ystep 1+ystep+xstep ...
           2+zstep 3+zstep 3+zstep+xstep 2+zstep+ystep];
    voxrow = bsxfun(@plus,(0:Xsz-1)'*xstep,vox)';
    voxplane = bsxfun(@plus,(0:Ysz-1)*ystep,voxrow(:));
    voxvol   = bsxfun(@plus,(0:Zsz-1)*zstep,voxplane(:));
    voxvol   = reshape(voxvol, 12, [])';
    I = kron((1:size(voxvol,1))',ones(1,12));
    nEdges  = 3*(Zsz+1)*(Ysz+1)*(Xsz+1);
    nVox    = Xsz*Ysz*Zsz;
    vox2edge = sparse(I,voxvol,ones(size(I)),nVox,nEdges);
    vol =      m.edge_length(voxvol(:,1)) ...
            .* m.edge_length(voxvol(:,2)) ...
            .* m.edge_length(voxvol(:,3));
    
    
%-------------------------------------------------------------------------%
% Show voxels
function show_voxels(rmdl,voxels)    
    mdl=rmdl;
    mdl.elems = voxels;
    mdl.boundary = mdl.elems;
    clf
    show_fem(mdl,[0 0 1]);

%-------------------------------------------------------------------------%
% Test indexing of faces on planes by showing random x, y and z plane
function test_faces(rmdl, faces, opt)
    mdl = rmdl;
    mdl.elems = faces;
    mdl.boundary = mdl.elems;
    img = mk_image(mdl,0);
    img.elem_data(opt.xplane + (randi(opt.Xsz+1)-1)*opt.xstep) = 1;
    img.elem_data(opt.yplane + (randi(opt.Ysz+1)-1)*opt.ystep) = 2;
    img.elem_data(opt.zplane + (randi(opt.Zsz+1)-1)*opt.zstep) = 3;
    show_fem(img,[0 0 1]);

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
    eidors_msg('models scaled by %g', scale,2);
    
%-------------------------------------------------------------------------%
% Parse option struct
 function opt = parse_opts(fmdl,rmdl, opt)

    opt.xvec = unique(rmdl.nodes(:,1));
    opt.yvec = unique(rmdl.nodes(:,2));
    opt.zvec = unique(rmdl.nodes(:,3));
    opt.Xsz  = numel(opt.xvec)-1;
    opt.Ysz  = numel(opt.yvec)-1;
    opt.Zsz  = numel(opt.zvec)-1;
    opt.xstep = 3;
    opt.ystep = opt.xstep*opt.Xsz+1;
    opt.zstep = opt.ystep*(opt.Ysz+1);
    xrow    = 1 + (0:opt.Ysz-1)*opt.ystep;
    opt.xplane  = bsxfun(@plus, (0:opt.Zsz-1)'*opt.zstep,xrow);
    yrow    = 2 + (0:opt.Xsz-1)*opt.xstep;
    opt.yplane  = bsxfun(@plus, (0:opt.Zsz-1)'*opt.zstep,yrow);
    zrow    = 3 + (0:opt.Xsz-1)*opt.xstep;
    opt.zplane  = bsxfun(@plus, (0:opt.Ysz-1)'*opt.ystep,zrow);
    opt.nVox = opt.Xsz*opt.Ysz*opt.Zsz;
    opt.nTet = size(fmdl.elems,1);
    
    if ~isfield(opt, 'tol_node2tet');
        opt.tol_node2tet = eps; % * max(rmdl_rng,fmdl_rng)^3;
    end
    if ~isfield(opt, 'tol_edge2edge')
        opt.tol_edge2edge = 2*sqrt(3)*eps(min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))));
    end
    if ~isfield(opt, 'tol_edge2tri')
        opt.tol_edge2tri = eps; %1e-10
    end
    if ~isfield(opt, 'save_memory')
       opt.save_memory = 0;
    end
    eidors_msg('@@@ node2tet  tolerance = %g', opt.tol_node2tet,2);
    eidors_msg('@@@ edge2edge tolerance = %g', opt.tol_edge2edge,2);
    eidors_msg('@@@ edge2tri  tolerance = %g', opt.tol_edge2tri,2);

%-------------------------------------------------------------------------%
% fprintf wrapper to use eidors log_level
function logmsg(varargin)
    if eidors_msg('log_level') >= 2
        fprintf(varargin{:});
    end


%-------------------------------------------------------------------------%
% fprintf wrapper to use eidors log_level
function progmsg(num, erase)
    if eidors_msg('log_level') < 2, return; end
    if nargin == 1, erase = true; end
    str = '%3d%%';
    if erase
        str = ['\b\b\b\b' str ];
    end
    fprintf(str, round(100*num));
    
%-------------------------------------------------------------------------%
% Perfom unit tests
function do_unit_test
%     do_small_test;
    do_realistic_test;
%     figure
%     do_case_tests;
%     do_edge2edge_timing_test;
    
    
function do_case_tests
    ll = eidors_msg('log_level');
    eidors_msg('log_level',1);
    vox = mk_grid_model([],0:1,0:1,0:1);
    tet.type = 'fwd_model';
    tet.elems = [1 2 3 4];

    X = 4; Y = 6;
    for i = 1:30
        tet.nodes = [0 0 0; 0 1 0; 1 0 0; 0 0 1];
        fprintf('%d\n',i);
        switch i
            case 1 % nothing in common
                txt = 'Nothing in common';
                tet.nodes = tet.nodes + 2;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                mk_grid_c2f(tet,vox);

            case 2 % common node
                tet.nodes = tet.nodes + 1;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp('Comm node tedge2rec',size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp('Comm node vedge2tri',size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp('Comm node edge2edge',size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp('Comm node tnode2vox',fnode2vox,[1;0;0;0],0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp('Comm node vnode2tet',nnz(rnode2tet),0,0);  % nothing
                mk_grid_c2f(tet,vox);

            case 3 % tet_edge v vox_node
                txt = 'tet_edge v vox_node';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) - 0.5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                res = zeros(8,1); res(1) = 1;
                unit_test_cmp(txt,rnode2tet,res,0);  % 1 point
                mk_grid_c2f(tet,vox);

            case 4 % tet_edge v vox_edge
                txt = 'tet_edge v vox_edge';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) - 0.5;
                tet.nodes(:,3)   = tet.nodes(:,3)   + 0.5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),1,0); % 1 point
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0); % nothing
                mk_grid_c2f(tet,vox);

            case 5 % tet_edge on vox_face
                txt = 'tet_edge on vox_face';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) - 0.3;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),1,0); % 1 point
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),2,0); % 2 points
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),1,0); % 1 point
                mk_grid_c2f(tet,vox);

            case 6 
                txt = 'vox_node on tet_face';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) - 0.3;
                tet.nodes(:,3)   = tet.nodes(:,3) -0.4;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),1,0); % 1 point
                mk_grid_c2f(tet,vox);

            case 7
                txt = 'tet_node on vox_face';
                tet.nodes(:,1) = tet.nodes(:,1) - 1;
                tet.nodes(:,2:3)  = tet.nodes(:,2:3) + .5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0); % nothing  
                mk_grid_c2f(tet,vox);


            case 8
                txt = 'tet_node on vox_edge';
                tet.nodes(:,1) = tet.nodes(:,1) - 1;
                tet.nodes(:,2)  = tet.nodes(:,2) + .5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0); % nothing   
                mk_grid_c2f(tet,vox);

            case 9
                txt = 'all nodes common';
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),4,0); % 4 points
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0); % nothing
                mk_grid_c2f(tet,vox);

            case 10
                txt = 'tet in vox';
                tet.nodes = .25 + .5*tet.nodes;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),4,0); % 4 points
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0); % nothing  
                mk_grid_c2f(tet,vox);
                c2f = mk_grid_c2f(tet,vox);
                max_vox = max(vox.nodes);
                min_vox = min(vox.nodes);
                edges = max_vox-min_vox;
                vox_vol = edges(:,1) .* edges(:,2) .* edges(:,3);
                tet_vol = get_elem_volume(tet);
                unit_test_cmp(txt,c2f, 1, 0);

            case 11
                txt = 'vox in tet';
                tet.nodes =  4*tet.nodes - .25;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0);  % nothing  
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),8,0); % 8 points  
                c2f = mk_grid_c2f(tet,vox);
                max_vox = max(vox.nodes);
                min_vox = min(vox.nodes);
                edges = max_vox-min_vox;
                vox_vol = edges(:,1) .* edges(:,2) .* edges(:,3);
                tet_vol = get_elem_volume(tet);
                unit_test_cmp(txt,c2f, vox_vol / tet_vol, 0);


            case 12 
                txt = 'tet_edge v. vox_face';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) - 0.3;
                tet.nodes(:,3)   = tet.nodes(:,3) + .4;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),2,0); % 2 points
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),2,0); % 2 points
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0); % nothing
                mk_grid_c2f(tet,vox);


            case 13
                txt = 'everything';
                tet.nodes = tet.nodes + 0.7;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),3,0); % 3 points
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),3,0); % 3 points
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),1,0); % 1 point
                mk_grid_c2f(tet,vox);

                %------- degenerate cases-----------%
            case 14 % common edge
                txt = 'DG1: Common edge';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) + 1;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,fnode2vox,[1;0;0;1],0); % 2 points
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0);  % nothing
                mk_grid_c2f(tet,vox);

            case 15 % common edge
                txt = 'DG2: Common edge';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) + 1;
                tet.nodes(:,3) = tet.nodes(:,3) + 0.5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),1,0);  % 1 point
                mk_grid_c2f(tet,vox);

             case 16 % common edge
                txt = 'DG3: Common edge';
                tet.nodes(:,1:2) = tet.nodes(:,1:2) + 1;
                z = tet.nodes(:,3) == 0;
                tet.nodes(z,3) = -.5;
                tet.nodes(~z,3) = 1.5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),2,0);  % 2 points
                mk_grid_c2f(tet,vox);

            case 17 % edge on face
                txt = 'DG4: edge on face';
                tet.nodes = [0 0 0; 1 .5 0; 0 1 0; 1 .5 1];
                tet.nodes(:,1) = tet.nodes(:,1) - 1;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),2,0); % 2 points
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0);  % nothing
                mk_grid_c2f(tet,vox);

            case 18 % edge2edge only
                txt = 'DG5: edge2edge only';
                tet.nodes = [0 0 0; 1 .5 0; 0 1 0; 1 .5 1];
                tet.nodes(:,1) = tet.nodes(:,1) - 1;
                z = tet.nodes(:,3) == 0;
                tet.nodes(z,3) = -.5;
                tet.nodes(~z,3) = 1.5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),2,0); % 2 points
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),0,0); % nothing
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0);  % nothing
                mk_grid_c2f(tet,vox);

            case 19 % edge on face
                txt = 'DG6: edge on face';
                tet.nodes = [0 0 0; 1 .5 0; 0 1 0; 1 .5 1];
                tet.nodes(:,1) = tet.nodes(:,1) - 1;
                tet.nodes(:,3) = tet.nodes(:,3) - .5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),1,0); % 1 point
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 points
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0);  % nothing
                mk_grid_c2f(tet,vox);

            case 20 %face on face
                txt = 'DG7: face on face';
                tet.nodes(:,1) = tet.nodes(:,1) + 1;
                tet.nodes(:,2) = tet.nodes(:,2) + 0.5;
                tet.nodes(:,3) = tet.nodes(:,3) + 0.5;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(unique(intpts,'rows'),1),2,0); % 2 unique points (these edges are counted 3 times, but vox2edge takes care of this)
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),1,0);  % 1 point
                mk_grid_c2f(tet,vox);

           case 21 %face on face
                txt = 'DG8: face on face';
                tet.nodes(:,1) = tet.nodes(:,1) + 1;
                tet.nodes(:,2) = tet.nodes(:,2) + 0.4;
                tet.nodes(:,3) = tet.nodes(:,3) - 0.6;
                subplot(X,Y,i), show_test(vox,tet);
                [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(intpts,1),0,0); % nothing
                [intpts, face2edge, face2intpt, edge2intpt] = test_tedge2vedge_intersections(vox,tet);
                unit_test_cmp(txt,size(unique(intpts,'rows'),1),2,0); % 2 unique points (these edges are counted 3 times, but vox2edge takes care of this)
                [fnode2vox] = test_tnode_in_vox(vox,tet);
                unit_test_cmp(txt,nnz(fnode2vox),1,0); % 1 point
                rnode2tet = test_vnode_in_tet(vox,tet);
                unit_test_cmp(txt,nnz(rnode2tet),0,0);  % nothing
                mk_grid_c2f(tet,vox);

            otherwise
                break;
        end
    end
    eidors_msg('log_level',ll);

function do_small_test
   fmdl= ng_mk_cyl_models([2,2,.1],[16,1],[.1,0,.025]);
   xvec = [-1.5:1:1.5];
   yvec = [-1.6:.8:1.6];
   zvec = -1:1:2;
   rmdl = mk_grid_model([],xvec,yvec,zvec);
%    profile clear
%    profile on
   eidors_cache off mk_grid_c2f
   tic
   c2f_a = mk_grid_c2f(fmdl, rmdl);
   toc
%    eidors_cache on mk_grid_c2f
   tic
   opt.save_memory = 1;
   c2f_b = mk_grid_c2f(fmdl, rmdl, opt);
   toc
%    size(c2f_a)
   assert(any(c2f_a(:) ~= c2f_b(:)) == 0);
%    profile off
%    p = profile('info');
%    profview(0,p);
%    fname = sprintf('mk_grid_c2f_profile_%s_%s_%s.mat',...
%        computer('arch'),version('-release'),datestr(now,'yyyymmdd-HHMMSS'));
%    tmp = tempname;
%    diary(tmp);
%    ver('Matlab');
%    diary off
%    str = fileread(tmp);
%    delete(tmp);
%    save(fname,'p','str');
%    eidors_msg('Profile file saved under %s',[cd filesep fname],0);
%    eidors_msg('Total time: %f s', max(cell2mat({p.FunctionTable.TotalTime})),0);
   
   
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
fprintf('Analytic: t=%f s\n',t);

tic
c2f_n = mk_coarse_fine_mapping(fmdl,rmdl);
t = toc;
fprintf('Approximate: t=%f s\n',t);


tetvol = get_elem_volume(fmdl);
opt = parse_opts(fmdl,rmdl);
m = prepare_vox_mdl(rmdl,opt);


f2c_a = bsxfun(@times, c2f_a, tetvol);
f2c_a = bsxfun(@rdivide,f2c_a', m.volume); 
img = mk_image(rmdl,0);
img.elem_data = f2c_a*ones(size(fmdl.elems,1),1);
subplot(132)
show_fem(img);

f2c_n = bsxfun(@times, c2f_n, tetvol);
f2c_n = bsxfun(@rdivide,f2c_n', m.volume); 
img = mk_image(rmdl,0);
img.elem_data = f2c_n*ones(size(fmdl.elems,1),1);
subplot(133)
show_fem(img);
subplot(131);
h = show_fem(fmdl);
set(h,'LineWidth',0.1)
hold on
h = show_fem(rmdl);
set(h,'EdgeColor','b','LineWidth',2);
hold off

function do_edge2edge_timing_test
    fmdl= ng_mk_cyl_models([2,2],[8,1],[.15,0,.05]);
    xvec = linspace(-2,2,33);
    yvec = linspace(-2,2,33);
    zvec = 0:.5:2;
    rmdl = mk_grid_model([],xvec,yvec,zvec);
    test_tedge2vedge_intersections(rmdl,fmdl);
    
    
function rnode2tet = test_vnode_in_tet(rmdl,fmdl)
    opt  = parse_opts(fmdl,rmdl);
    fmdl = prepare_fmdl(fmdl);
    rnode2tet = get_nodes_in_tets(fmdl,rmdl, opt);

function [fnode2vox] = test_tnode_in_vox(rmdl,fmdl)
    fmdl = prepare_fmdl(fmdl);
    [fnode2vox] = get_nodes_in_voxels(fmdl,rmdl);
    
function [intpts, face2edge, face2intpt, edge2intpt] = test_rec_tedge_intersections(rmdl,fmdl)
    opt = parse_opts(fmdl,rmdl);
    m = prepare_vox_mdl(rmdl,opt);
    fmdl = prepare_fmdl(fmdl);
    [intpts, face2edge, face2intpt, edge2intpt] = get_voxel_intersection_points(fmdl,m.faces,opt);
    
    
function [intpts, face2edge, face2intpt, edge2intpt] = test_tet_vedge_intersections(rmdl,fmdl)
    opt = parse_opts(fmdl,rmdl);
    m = prepare_vox_mdl(rmdl,opt);
    fmdl = prepare_fmdl(fmdl);
    [intpts, face2edge, face2intpt, edge2intpt] = get_tet_intersection_points(fmdl,m,opt);
    
function [intpts, tedge2vedge, tedge2intpt, vedge2intpt] = test_tedge2vedge_intersections(rmdl,fmdl)
    opt = parse_opts(fmdl,rmdl);
    m = prepare_vox_mdl(rmdl,opt);
    fmdl = prepare_fmdl(fmdl);
    [intpts, tedge2vedge, tedge2intpt, vedge2intpt] = find_edge2edge_intersections(fmdl.edges,fmdl.nodes,m.edges,rmdl.nodes, opt.tol_node2tet);

function show_test(vox,tet)
    show_fem(vox);
    hold on
    h = show_fem(tet);
    set(h, 'EdgeColor','b');
    hold off
    axis auto
