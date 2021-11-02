function point2tet = point_in_tet(fmdl,points, epsilon)
%POINT_IN_TET test for points contained in elements

% (C) 2015 Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

copt.str = 'point_in_tet';
copt.cache_obj = {fmdl.nodes, fmdl.elems, points, epsilon};

point2tet = eidors_cache(@do_point_in_tet,{fmdl, points, epsilon}, copt);

end

function point2tet = do_point_in_tet(fmdl, points, epsilon)
   progress_msg('Point_in_tet');
   [A,b] = tet_to_inequal(fmdl.nodes,fmdl.elems);
   progress_msg(.01);
   % This is split to decrease the memory footprint
   point2tet = (bsxfun(@minus, A(1:4:end,:)*points',b(1:4:end)) <= epsilon)';
   progress_msg(.21);
   for i = 2:4
      point2tet = point2tet & (bsxfun(@minus, A(i:4:end,:)*points',b(i:4:end)) <= epsilon)';
      progress_msg(.21 + (i-1)*.23);
   end

   % exclude coinciding nodes
   ex= bsxfun(@eq,points(:,1),fmdl.nodes(:,1)') & ...
       bsxfun(@eq,points(:,2),fmdl.nodes(:,2)') & ...
       bsxfun(@eq,points(:,3),fmdl.nodes(:,3)');
   progress_msg(.94);
   point2tet(any(ex,2),:) = 0;
   point2tet = sparse(point2tet);
   progress_msg(Inf);
end