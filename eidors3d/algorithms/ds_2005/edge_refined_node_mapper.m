function [index_vtx] = vex(vtx_coarse,vtx_dense,simp_coarse,simp_dense);

% This function maps a dense mesh verticies array onto a more coarse mesh verticies
% array.  The closest vertex on the dense mesh to the objective vertex on
% the coarse mesh is found.
%
% vtx_coarse  = vertices array from the coarse mesh
% vtx_dense   = vertices array from the dense mesh
% simp_coarse = simplicies array from the coarse mesh
% simp_dense  = simplicies array from the dense mesh
% index_vtx   = the index array mapping each coarse mesh (from netgen) vtx onto a dense mesh (from netgen) vtx

% Array pre-allocation

dist=zeros(size(vtx_dense,1),size(vtx_coarse,1));
index=zeros(size(vtx_coarse,1),2);


% Down to business ...

h = waitbar(0,'Calculating Verticies Map');

for ic=1:size(vtx_coarse,1);   % for all coarse verticies
    
    waitbar(ic/size(vtx_coarse,1))


    for id=1:size(vtx_dense,1);   % for all dense verticies
    
        dx=vtx_dense(id,1)-vtx_coarse(ic,1);   % find the x co-ord difference
        dy=vtx_dense(id,2)-vtx_coarse(ic,2);   % find the y co-ord difference
        dz=vtx_dense(id,3)-vtx_coarse(ic,3);   % find the z co-ord difference
    
        dist(id,ic)=sqrt((dx^2)+(dy^2)+(dz^2));   % distance between points for each dense vertex and the ic'th coarse vertex
       
        id=id+1;
        
    end 
    
    [m,index(ic,1)]=min(dist(:,ic));   % index out the minimum distance from the dense mesh to the ic'th vertex
    
    index(ic,2)=m;   % write the actual minimum distance (as a quality control procedure)
    
    ic=ic+1;
    
end

close(h)

index_vtx=index;

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