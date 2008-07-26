function [srf] = find_boundary(simp);
% [srf] = find_boundary(simp);
%
%Caclulates the boundary faces of a given 3D volume.
%Usefull in electrode assignment.
%
%srf  =  array of elements on each boundary simplex
%        boundary simplices are of 1 lower dimention than simp
%simp = The simplices matrix

% $Id$

wew = size(simp,2) - 1;

if wew==3 || wew==2
   srf= find_2or3d_boundary(simp,wew);
else
   eidors_msg('find_boundary: WARNING: not 2D or 3D simplices',1);
   srf=[]; return;
end

% sort surfaces. If there is more than one, its not on the boundary
function srf= find_2or3d_boundary(simp,dim);
   % convert to integer to make sort faster
   simp = uint32( simp ); % max 4b elements - enough for a while
   localface = nchoosek(1:dim+1,dim);
   srf_local= simp(:,localface');
   srf_local= reshape( srf_local', dim, []);
   srf_local= sort(srf_local)';
   sort_srl = sortrows( srf_local );
   same_srl = find( all( sort_srl(1:end-1,:) == sort_srl(2:end,:), 2) );
   diff_srl = ones(size(srf_local,1),1);
   diff_srl(same_srl) = 0;
   diff_srl(same_srl+1) = 0;

   srf= sort_srl( find(diff_srl),: );

