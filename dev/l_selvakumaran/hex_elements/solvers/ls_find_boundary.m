function [srf, idx] = ls_find_boundary(simp);
%Finds boundary of mesh with quad/hex elements

wew = size(simp,2);

if  wew==4; 
    dim=2;
   [srf,idx]= find_2d_boundary(simp,dim);
elseif wew==8;
    dim=3;
   [srf,idx]= find_3d_boundary(simp,dim);
else
   eidors_msg('find_boundary: WARNING: not 2D or 3D simplices',1);
   srf=[]; return;
end

function [srf,idx]= find_2d_boundary(simp,dim);
   if size(simp,1) < 4e9; % max of uint32
      % convert to integer to make sort faster
      simp = uint32( simp );
   end
   %localface = nchoosek(1:dim+1,dim)
   localface=[1 2;2 3;3 4;4 1];
   srf_local= simp(:,localface');
   srf_local= reshape( srf_local', 2, []); % D x 3E
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
   idx= ceil(idx/(4));
   
   function [srf,idx]= find_3d_boundary(simp,dim);
   if size(simp,1) < 4e9; % max of uint32
      % convert to integer to make sort faster
      simp = uint32( simp );
   end
   %localface = nchoosek(1:dim+1,dim)
   localface=[1 2 3 4;1 2 6 5;1 5 8 4;2 6 7 3; 3 7 8 4; 5 6 7 8];
   srf_local= simp(:,localface');
   srf_local= reshape( srf_local', 4, []); % D x 3E
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
   idx= ceil(idx/(6));