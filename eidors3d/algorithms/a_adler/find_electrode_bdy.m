function [bdy_idx, bdy_area] = find_electrode_bdy( bdy, vtx, elec_nodes)
% FIND_ELECTRODE_BDY: find the boundary index are area for electrode
% Input:
%   bdy => boundary (from fwd_model.boundary) bdy simplices x nodes index 
%   vtx => node pts (from fwd_model.nodes)
%   elec_nodes => index of nodes in the electrode
% Output:
%   bdy_idx  => vector of boundary simplices in this electrode
%   bdy_area => boundary area of each simplex in bdy_idx

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: find_electrode_bdy.m,v 1.1 2008-05-14 21:33:58 aadler Exp $

bdy_idx = find_bdy_idx( bdy, elec_nodes);

l_bdy_idx= length(bdy_idx);
bdy_area= zeros(1, l_bdy_idx);

for i=1:l_bdy_idx
   bdy_nodes   = bdy(bdy_idx(i),:);
   bdy_area(i) = tria_area( vtx(bdy_nodes,:) );
end


function ffb = find_bdy_idx( bdy, elec_nodes);
   bdy_els = zeros(size(bdy,1),1);
   for nd= unique(elec_nodes);
      bdy_els = bdy_els + any(bdy==nd,2);
   end
   ffb = find(bdy_els == size(bdy,2));

% bdy points is [x1,y1,z1;x2,y2,z2; etc]
function area= tria_area( bdy_pts ); 
   vectors= diff(bdy_pts); 
   if size(vectors,1)==2
      vectors= cross(vectors(1,:),vectors(2,:));
   end
   d= size(bdy_pts,1);
   area= sqrt( sum(vectors.^2) )/( d-1 );

