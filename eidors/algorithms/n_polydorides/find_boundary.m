function [srf, idx] = find_boundary(simp);
% [srf] = find_boundary(simp);
%
%Caclulates the boundary faces of a given 3D volume.
%Usefull in electrode assignment.
%
%srf  =  array of elements on each boundary simplex
%        boundary simplices are of 1 lower dimention than simp
%idx  =  index of simplex to which each boundary belongs
%simp = The simplices matrix

% $Id$

wew = size(simp,2) - 1;

if wew==3 || wew==2
   [srf,idx]= find_2or3d_boundary(simp,wew);
else
   eidors_msg('find_boundary: WARNING: not 2D or 3D simplices',1);
   srf=[]; return;
end

% sort surfaces. If there is more than one, its not on the boundary
function [srf,idx]= find_2or3d_boundary(simp,dim);
   if size(simp,1) < 4e9 % max of uint32
      % convert to integer to make sort faster
      simp = uint32( simp );
   end
   localface = nchoosek(1:dim+1,dim);
   srf_local= simp(:,localface');
   srf_local= reshape( srf_local', dim, []); % D x 3E
   srf_local= sort(srf_local)'; % Sort each row
   [sort_srl,sort_idx] = sortrows( srf_local );

   % Fine the ones that are the same
   first_ones =  sort_srl(1:end-1,:);
   next_ones  =  sort_srl(2:end,:);
   same_srl = find( all( first_ones == next_ones, 2) );

   % Assume they're all different. then find the same ones
   diff_srl = logical(ones(size(srf_local,1),1));
   diff_srl(same_srl) = 0;
   diff_srl(same_srl+1) = 0;

   srf= sort_srl( diff_srl,: );
   idx= sort_idx( diff_srl);
   idx= ceil(idx/(dim+1));


