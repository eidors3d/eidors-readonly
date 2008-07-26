function [bdy_idx, bdy_area] = find_electrode_bdy( bdy, vtx, elec_nodes)
% FIND_ELECTRODE_BDY: find the boundary index are area for electrode
% Input:
%   bdy => boundary (from fwd_model.boundary) bdy simplices x nodes index 
%   vtx => node pts (from fwd_model.nodes)
%   elec_nodes => index of nodes in the electrode
% Output:
%   bdy_idx  => vector of boundary simplices in this electrode
%   bdy_area => boundary area of each simplex in bdy_idx
%  if the nodes in the electrode are points, then
%   bdy_area is the area corresponding to these point electrodes

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

[bdy_idx, point] = find_bdy_idx( bdy, elec_nodes);
if nargout==1; return;end

l_bdy_idx= length(bdy_idx);
l_point= length(point);

if l_bdy_idx > 0 && l_point ==0
   bdy_area= zeros(1, l_bdy_idx);

   for i=1:l_bdy_idx
      bdy_nodes   = bdy(bdy_idx(i),:);
      bdy_area(i) = tria_area( vtx(bdy_nodes,:) );
   end
elseif l_bdy_idx == 0 && l_point >0
   dims = size(bdy,2); 
   l_point= length(point);
   bdy_area = zeros(l_point,1);
   for i=1:l_point
     ff= find( any(bdy== point(i),2) );
     this_area= 0;
     for ffp=ff(:)'
        xyz= vtx( bdy(ffp,:),:);
        this_area= this_area + tria_area( xyz );
     end
     bdy_area(i)= bdy_area(i) + this_area/dims;
   end
else
   error('can`t model this electrode, with %d CEM and %d point', ...
            l_bdy_idx, l_point )
end
   


% Find which boundary nodes are part of a 
% complete simplex. These nodes can be made
% into complete electrode model elements
%
% Any nodes left will be treated as point electrodes
function [ffb,unused] = find_bdy_idx( bdy, elec_nodes);
   bdy_els = zeros(size(bdy,1),1);
   elec_nodes = unique(elec_nodes);
   for nd= elec_nodes(:)'
      bdy_els = bdy_els + any(bdy==nd,2);
   end
   ffb = find(bdy_els == size(bdy,2));

%  find if all nodes are used
   used_nodes= reshape( bdy(ffb,:), 1,[]);
   unused=     setdiff( elec_nodes, used_nodes);
   

% bdy points is [x1,y1,z1;x2,y2,z2; etc]
function area= tria_area( bdy_pts ); 
   vectors= diff(bdy_pts); 
   if size(vectors,1)==2
      vectors= cross(vectors(1,:),vectors(2,:));
   end
   d= size(bdy_pts,1);
   area= sqrt( sum(vectors.^2) )/( d-1 );

