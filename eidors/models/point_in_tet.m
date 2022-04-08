function point2tet = point_in_tet(fmdl,points, epsilon, exclude_nodes)
%POINT_IN_TET test for points contained in elements
%
%
% (C) 2015 Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

if nargin < 4
    exclude_nodes = false;
end

copt.fstr = 'point_in_tet';
copt.cache_obj = {fmdl.nodes, fmdl.elems, points, epsilon, exclude_nodes};

point2tet = eidors_cache(@do_point_in_tet,{fmdl, points, epsilon, exclude_nodes}, copt);

end

function point2tet = do_point_in_tet(fmdl, points, epsilon, exclude_nodes)
   ver = eidors_obj('interpreter_version');
   if ~ver.isoctave 
       warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
       TR = triangulation(fmdl.elems, fmdl.nodes);
       warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
       ID = pointLocation(TR, points);
       idx = ~isnan(ID);
       point2tet = builtin('sparse',find(idx),ID(idx),1,size(points,1),size(fmdl.elems,1));
   else
       progress_msg('Tet inequalities');
       [A,b] = tet_to_inequal(fmdl.nodes,fmdl.elems);
       progress_msg(Inf);
       
       %split computation into managable chunks to fit in memory
       mem_req = size(points,1) * size(fmdl.elems,1) * 8; % in bytes
       mem_use = 2*(1024^3); % 2 GiB
       n_chunks = ceil(mem_req / mem_use);
       chunk_sz = ceil(size(points,1) / n_chunks);
       point2tet = logical(builtin('sparse',0, size(fmdl.elems,1)));
       chunk_end = 0;
       progress_msg('Point in tet',0,n_chunks);
       for c = 1:n_chunks
           progress_msg(c,n_chunks)
           chunk_start = chunk_end+1;
           chunk_end = min(chunk_start + chunk_sz, size(points,1));
           idx = chunk_start:chunk_end;
           if 1
               p2t = (bsxfun(@minus, A(1:4:end,:)*points(idx,:)',b(1:4:end)) <= epsilon)';
               %        progress_msg(.21);
               for i = 2:4
                   % that was actually slower
                   %good = find(any(p2t,2));
                   %p2t(good,:) = p2t(good,:) & (bsxfun(@minus, A(i:4:end,:)*points(idx(good),:)',b(i:4:end)) <= epsilon)';
                   p2t = p2t & (bsxfun(@minus, A(i:4:end,:)*points(idx,:)',b(i:4:end)) <= epsilon)';
                   %           progress_msg(.21 + (i-1)*.23);
               end
               point2tet = [point2tet; builtin('sparse',p2t)];
           else
               % slower...
               p2t = (bsxfun(@minus, A*points(idx,:)',b) <= epsilon)';
               point2tet = [point2tet; builtin('sparse',reshape(all(reshape(p2t',4,[])),[],length(idx))')];
           end
       end
       progress_msg(Inf);
   end
       % exclude coinciding nodes
       %    ex= bsxfun(@eq,points(:,1),fmdl.nodes(:,1)') & ...
       %        bsxfun(@eq,points(:,2),fmdl.nodes(:,2)') & ...
       %        bsxfun(@eq,points(:,3),fmdl.nodes(:,3)');
       %    progress_msg(.94);
       %    point2tet(any(ex,2),:) = 0;
       
   if exclude_nodes    
       ex = ismember(points, fmdl.nodes, 'rows');
       point2tet(ex,:) = 0;
   end
       

end