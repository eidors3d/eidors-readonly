function [bdy_idx, bdy_area] = ls_find_electrode_bdy( bdy, vtx, elec_nodes)
% FIND_ELECTRODE_BDY: find the boundary index area for electrode when the
% mesh contains hex/quad elements
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

if l_bdy_idx > 0 && l_point ==0;
   bdy_area= zeros(1, l_bdy_idx);

   for i=1:l_bdy_idx;
      bdy_nodes   = bdy(bdy_idx(i),:);
      bdy_area(i) = tria_area( vtx(bdy_nodes,:) );
   end
elseif l_bdy_idx == 0 && l_point >0;
   dims = size(bdy,2); 
   l_point= length(point);
   bdy_area = zeros(l_point,1);
   for i=1:l_point;
     ff= find( any(bdy== point(i),2) );
     this_area= 0;
     for ffp=ff(:)';
        xyz= vtx( bdy(ffp,:),:);
        this_area= this_area + tria_area( xyz );
     end
     bdy_area(i)= bdy_area(i) + this_area/dims;
   end
else
   error('Can''t model this electrode. It has %d CEM and %d point nodes.', ...
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
   for nd= elec_nodes(:)';
      bdy_els = bdy_els + any(bdy==nd,2);
   end
   ffb = find(bdy_els == size(bdy,2));

%  find if all nodes are used
   used_nodes= reshape( bdy(ffb,:), 1,[]);
   unused=     setdiff( elec_nodes, used_nodes);
   

% bdy points is [x1,y1,z1;x2,y2,z2; etc]
function area= tria_area( bdy_pts ); 
   if size(bdy_pts,1)==2;
       area=sqrt( sum(diff(bdy_pts).^2) ); 
   elseif size(bdy_pts,1)==4;
       tri_1=bdy_pts(1:3,:);
       tri_2=[bdy_pts(3,:);bdy_pts(4,:);bdy_pts(1,:)];
       vectors_1=diff(tri_1);
       cross_1= cross(vectors_1(1,:),vectors_1(2,:));
       area_1= sqrt( sum(cross_1.^2) )/2;
       vectors_2=diff(tri_2);
       cross_2= cross(vectors_2(1,:),vectors_2(2,:));
       area_2= sqrt( sum(cross_2.^2) )/2;
       area=area_1+area_2;
   end
   
   

