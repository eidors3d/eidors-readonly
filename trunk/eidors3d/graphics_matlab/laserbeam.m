function [sel] = laserbeam(vtx,srf,cnts)
% function [sel] = laserbeam(vtx,srf,cnts)
%
% Auxiliary ploting function, to assist in assigning surface patches
% as electrodes. 
%
%
%
%vtx  = The vertices matrix
%srf  = The boundary faces 
%cnts = The coordinates of the centers of the surfaces


p = get(gca,'CurrentPoint');

pnear = p(1,:);

  
   x1 = pnear(1)*ones(size(cnts,1),1);     
   y1 = pnear(2)*ones(size(cnts,1),1);     
   z1 = pnear(3)*ones(size(cnts,1),1);
   x2 = cnts(:,1); 
   y2 = cnts(:,2); 
   z2 = cnts(:,3);
   
[dista] = db23d(x1,y1,z1,x2,y2,z2);


[mv,sel] = min(dista);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
