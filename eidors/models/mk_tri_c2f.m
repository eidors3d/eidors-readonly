function c2f = mk_tri_c2f(fmdl,rmdl,opt)
%MK_TRI_C2F

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


copt.fstr = 'mk_tri_c2f';

c2f(fmdl_idx,rmdl_idx) = eidors_cache(@do_mk_tri_c2f,{fmdl,rmdl,opt},copt);
end

function c2f = do_mk_tri_c2f(fmdl,rmdl,opt)
   DEBUG = eidors_debug('query','mk_tri_c2f');
   
   c2f = sparse(0,0);
   progress_msg('Prepare fine model...');
   fmdl = prepare_tri_mdl(fmdl);
   progress_msg(Inf);
   
   progress_msg('Prepare course model...');
   rmdl = prepare_tri_mdl(rmdl);
   progress_msg(Inf);
   
   progress_msg('Find edge2edge intersections...');
   [intpts,fedge2redge,fedge2intpts,redge2intpts] = ...
      find_edge2edge_intersections( fmdl.edges,fmdl.nodes,...
                                    rmdl.edges,rmdl.nodes, opt.tol_edge2edge);
   progress_msg(sprintf('Found %d',size(intpts,1)),Inf);

   
   progress_msg('Find c_nodes in f_elems...');
   rnode2ftri = point_in_triangle(rmdl.nodes, fmdl.elems, fmdl.nodes, opt.tol_node2tri);
   progress_msg(sprintf('Found %d', nnz(rnode2ftri)), Inf);
   
   progress_msg('Find c_elems in f_elems...')
   rtri_in_ftri = (double(rmdl.node2elem') * rnode2ftri) == 3;
   progress_msg(sprintf('Found %d',nnz(rtri_in_ftri)), Inf);
   
   progress_msg('Find f_nodes in c_elems...');
   fnode2rtri = point_in_triangle(fmdl.nodes, rmdl.elems, rmdl.nodes, opt.tol_node2tri);
   progress_msg(sprintf('Found %d', nnz(fnode2rtri)), Inf);
   
   progress_msg('Find f_elems in c_elems...')
   ftri_in_rtri = (double(fmdl.node2elem') * fnode2rtri) == 3;
   progress_msg(sprintf('Found %d',nnz(ftri_in_rtri)), Inf);

   progress_msg('Find total intersections...');
   rtri2ftri = double(rmdl.edge2elem') * fedge2redge' * fmdl.edge2elem;
   % exclude inclusion (dealt with separately)
   rtri2ftri = rtri2ftri & ~rtri_in_ftri & ~ftri_in_rtri'; 
   progress_msg(sprintf('Found %d',nnz(rtri2ftri)), Inf);
   
   progress_msg('Calculate intersection volumes...');
   % sparse logical multiplication doesn't exist
   rtri2intpt = logical(rmdl.edge2elem'*redge2intpts)';
   ftri2intpt = logical(fmdl.edge2elem'*fedge2intpts)';
   
   rtri_todo = find(sum(rtri2ftri,2)>0);
   C = []; F = []; V = [];

   id = 0; N = length(rtri_todo);
   mint = ceil(N/100);
   for v = rtri_todo'
      id = id+1;
      if mod(id,mint)==0, progress_msg(id/N); end
      tri_todo = find(rtri2ftri(v,:));

      common_intpts = bsxfun(@and,rtri2intpt(:,v), ftri2intpt(:,tri_todo));
      
      f_nodes     = bsxfun(@and,fnode2rtri(:,v), fmdl.node2elem(:,tri_todo));
      r_nodes     = bsxfun(@and,rnode2ftri(:,tri_todo), rmdl.node2elem(:,v));
      
      C = [C; v*ones(numel(tri_todo),1)];
      F = [F; tri_todo'];
      last_v = numel(V);
      V = [V; zeros(numel(tri_todo),1)]; % pre-allocate
      
      for t = 1:numel(tri_todo)
         pts = [ intpts(common_intpts(:,t),:);
                  fmdl.nodes(f_nodes(:,t),:);
                  rmdl.nodes(r_nodes(:,t),:)];
         last_v = last_v + 1;
         if size(pts,1) < 3, continue, end 
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
            V(last_v) = max(V(last_v),0); % numerical issues may produce tiny negative volume
            V(last_v) = V(last_v) * scale^2; % undo scaling
         catch err
            ok = false;
            switch err.identifier
               case {'MATLAB:qhullmx:DegenerateData', 'MATLAB:qhullmx:UndefinedError'}
                  % border case point may be included multiple times.
                  % this is OK... though we may miss cases where more
                  % points should have been found but were not
                  u = uniquetol(pts*scale,eps,'ByRows',true,'DataScale', 1);
                  ok = ok | (size(u,1) < 3);
                  if ~ok
                     % test for colinearity
                     cp = bsxfun(@minus, u(2:end,:), u(1,:));
                     l = sqrt(sum(cp.^2,2));
                     cp = bsxfun(@rdivide, cp, l);
                     u = uniquetol(cp,eps,'ByRows',true,'DataScale',1);
                     ok = ok | size(u,1) == 1; % co-linear points
                  end
            end
            if ~ok
               if DEBUG || eidors_debug('query','mk_tet_c2f:convhulln');
                  tri.nodes = fmdl.nodes;
                  vox.nodes = rmdl.nodes;
                  tri.type = 'fwd_model';
                  vox.type = 'fwd_model';
                  vox.elems = rmdl.elems(v,:);
                  vox.boundary = vox.elems;
                  tri.elems = fmdl.elems(tri_todo(t),:);
                  clf
                  show_fem(vox)
                  hold on
                  h = show_fem(tri);
                  set(h,'EdgeColor','b')
                  pts = bsxfun(@plus,pts*scale,ctr);
                  plot(pts(:,1),pts(:,2),'o');
                  hold off
                  axis auto
                  keyboard
               else
                  fprintf('\n');
                  eidors_msg(['convhulln has thrown an error. ' ...
                     'Enable eidors_debug on mk_tri_c2f and re-run to see a debug plot'],0);
                  rethrow(err);
               end
            end
         end
         
      end
   end
   progress_msg(Inf);
   
   c2f = sparse(F,C,V,size(fmdl.elems,1),size(rmdl.elems,1));
   
   % add rtri contained in ftri
   try rmdl = rmfield(rmdl,'coarse2fine'); end % messes with volume
   c2f = c2f + bsxfun(@times, sparse(rtri_in_ftri), get_elem_volume(rmdl))';
   
   % normalize to fine volume
   vol = get_elem_volume(fmdl);
   c2f = bsxfun(@rdivide,c2f,vol);
   
   % add fine elems contained in coarse elems
   c2f = c2f + ftri_in_rtri;

end


function [pts,FE2CE,FE2pts,CE2pts] = find_edge2edge_intersections(FE,FN,CE,CN, epsilon)
   P1 = FN(FE(:,1),:);
   P2 = FN(FE(:,2),:);
   P3 = CN(CE(:,1),:);
   P4 = CN(CE(:,2),:);
   
   P21 = P2 - P1;
   P43 = P4 - P3;
   
   invden = ( bsxfun(@times, P21(:,1), P43(:,2)') - ...
              bsxfun(@times, P21(:,2), P43(:,1)')       ).^-1;
   P13_x = bsxfun(@minus,P1(:,1),P3(:,1)');
   P13_y = bsxfun(@minus,P1(:,2),P3(:,2)');
   s = ( bsxfun(@times,-P21(:,2), P13_x) + ...
         bsxfun(@times, P21(:,1), P13_y)) .* invden;
   t = ( bsxfun(@times,-P43(:,2)',P13_x) + ...
         bsxfun(@times, P43(:,1)',P13_y)) .* invden;

   FE2CE= s >= 0 & s <= 1 & t >= 0 & t <= 1;
   
   [fe, ce] = find(FE2CE);
   N_pts = size(fe,1);
   pts = zeros(N_pts,2);
   for i = 1:N_pts
      pts(i,:) = P1(fe(i),:) + t(fe(i),ce(i)) * P21(fe(i),:);
   end
   FE2CE = sparse(FE2CE);
   N_ce = size(CE,1);
   N_fe = size(FE,1);
   FE2pts = sparse(fe, 1:N_pts, ones(N_pts,1), N_fe, N_pts);
   CE2pts = sparse(ce, 1:N_pts, ones(N_pts,1), N_ce, N_pts);


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
% Remove obviously non-overlapping parts of the models
function [fmdl,rmdl,fmdl_idx,rmdl_idx] = crop_models(fmdl,rmdl)

% instead of the usual approach, we could calculate the convex hull and
% then use inpolygon... but that would require trusting inpolygon

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
   re_gt = any(reshape(all(reshape(r_gt(rmdl.elems',:),3,[])),[],2),2);
   re_lt = any(reshape(all(reshape(r_lt(rmdl.elems',:),3,[])),[],2),2);
   fe_gt = any(reshape(all(reshape(f_gt(fmdl.elems',:),3,[])),[],2),2);
   fe_lt = any(reshape(all(reshape(f_lt(fmdl.elems',:),3,[])),[],2),2);
   
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
   fmdl.elems = reshape(f_idx(f_n),[],3);
   
   rmdl.nodes = rmdl.nodes(r_used_nodes,:);
   fmdl.nodes = fmdl.nodes(f_used_nodes,:);
   
   % for the benefit of any (debug) plots later on
   rmdl = rmfield(rmdl,'boundary');
   fmdl = rmfield(fmdl,'boundary');
    
end

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
end

%-------------------------------------------------------------------------%
% Parse option struct
 function opt = parse_opts(fmdl,rmdl, opt)

    
    if ~isfield(opt, 'tol_node2tri');
        opt.tol_node2tri = eps; % * max(rmdl_rng,fmdl_rng)^3;
    end
    if ~isfield(opt, 'tol_edge2edge')
        opt.tol_edge2edge = 2*sqrt(3)*eps(min(max(abs(fmdl.nodes(:))),max(abs(rmdl.nodes(:)))));
    end
%     if ~isfield(opt, 'tol_node2tri')
%         opt.tol_edge2tri = eps; %1e-10
%     end
%     if ~isfield(opt, 'save_memory')
%        opt.save_memory = 0;
%     end
    eidors_msg('@@ node2tri  tolerance = %g', opt.tol_node2tri,2);
    eidors_msg('@@ edge2edge tolerance = %g', opt.tol_edge2edge,2);
%     eidors_msg('@@ edge2tri  tolerance = %g', opt.tol_edge2tri,2);
 end   

function do_unit_test
   do_small_test
end

function do_small_test
   imdl = mk_common_model('a2c',16);
   rmdl = imdl.fwd_model;
%    rmdl = ng_mk_cyl_models([0 .5],[],[]);
%    rmdl.nodes(:,1) = rmdl.nodes(:,1) + .5;
%    show_fem(rmdl,[0 0 1])
%    hold all
   imdl = mk_common_model('d2c',16);
   fmdl = imdl.fwd_model;
%    fmdl = ng_mk_cyl_models([0 .5 .1],[],[]);
%    h = show_fem(fmdl);
%    set(h,'edgecolor','b','facecolor','none');
%    hold off
%    axis auto
   c2f = mk_tri_c2f(fmdl,rmdl);
   
   
   clf
   img1 = mk_image(fmdl,0);
   img1.elem_data = sum(c2f,2);
   img1.calc_colous.ref_level = .5;
   img1.calc_colours.clim = .5;
   
   subplot(121)
   show_fem(img1);
   img2 = mk_image(rmdl,0);
   img2.elem_data = (c2f' * get_elem_volume(fmdl)) ./ get_elem_volume(rmdl);
   img2.calc_colous.ref_level = .5;
   img2.calc_colours.clim = .55;
   subplot(122)
   show_fem(img2);
   
   unit_test_cmp('Check C2F size', size(c2f),[length(fmdl.elems), length(rmdl.elems)]);
   unit_test_cmp('Check C2F max==1', max(c2f(:)), 1);
   unit_test_cmp('Check C2F min==0', min(c2f(:)), 0);
   
   f2c = bsxfun(@rdivide, bsxfun(@times, c2f, get_elem_volume(fmdl))', get_elem_volume(rmdl));
   unit_test_cmp('Check F2C max==1', max(sum(f2c,2)), 1, 1e-15);
   unit_test_cmp('Check F2C min==0', min(f2c(:)), 0);
end
