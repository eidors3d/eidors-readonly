function [index_vtx] = edge_refined_node_mapper(mdl_coarse, mdl_dense);
% EDGE_REFINED_NODE_MAPPER:
%      maps a dense mesh verticies array onto a more coarse mesh verticies
%      array.  The closest vertex on the dense mesh to the objective vertex on
%      the coarse mesh is found.
%
% Usage:
%  [index_vtx] = edge_refined_node_mapper(mdl_coarse, mdl_dense);
%
% mdl_coarse  = fwd_model of coarse mesh
% mdl_dense   = fwd_model of dense mesh

warning('EIDORS:deprecated','EDGE_REFINED_NODE_MAPPER is deprecated as of 06-Jun-2012. ');

index_vtx = eidors_obj('get-cache', mdl_dense, 'index_vtx', mdl_coarse);
if ~isempty(index_vtx)
   eidors_msg('edge_refined_node_mapper: using cached value', 2);
   return
end

vtx_dense  = mdl_dense.nodes;
simp_dense = mdl_dense.elems;
vtx_coarse  = mdl_coarse.nodes;
simp_coarse = mdl_coarse.elems;

index=zeros(size(vtx_coarse,1),2);


% Down to business ...

h = waitbar(0,'Calculating Verticies Map');

for ic=1:size(vtx_coarse,1);   % for all coarse verticies

   waitbar(ic/size(vtx_coarse,1))

   dx=vtx_dense(:,1)-vtx_coarse(ic,1);   % find the x co-ord difference
   dy=vtx_dense(:,2)-vtx_coarse(ic,2);   % find the y co-ord difference
   dz=vtx_dense(:,3)-vtx_coarse(ic,3);   % find the z co-ord difference

   % distance between points for each dense vertex and the ic'th coarse vertex
   dist=sqrt((dx.^2)+(dy.^2)+(dz.^2));

   [m,index(ic,1)]=min(dist);   % index out the minimum distance from the dense mesh to the ic'th vertex

   index(ic,2)=m;   % write the actual minimum distance (as a quality control procedure)

end

close(h)

index_vtx=index;

% Cache the restult - it depends on both dense and coarse mdl
eidors_obj('set-cache', mdl_dense, 'index_vtx', index_vtx, mdl_coarse);
eidors_msg('edge_refined_node_mapper: setting cached value', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) D.R Stephenson 2004
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version XXX
% MATLAB Version 6.5.0.180913a (R13)
% MATLAB License Number: 1560
% Operating System: Microsoft Windows XP Version 5.1 (Build 2600: Service Pack 1)
% Java VM Version: Java 1.3.1_01 with Sun Microsystems Inc. Java HotSpot(TM) Client VM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
