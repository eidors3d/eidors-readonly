function potplot(vtx,srf,V);
%function potplot(vtx,srf,V);
%
%Function that animates the forward solution in 3D for all existing 
%current patterns.
%
%
%
%vtx = The vertices matrix
%srf = The boundary surfaces
%V   = The forward solution

figure;

for i=1:size(V,2);
   
Vp = V(:,i); 
   
trisurf(srf,vtx(:,1),vtx(:,2),vtx(:,3),Vp(1:size(vtx,1)));
shading interp;
daspect([1 1 1]);

view(3);
camzoom(1);
grid off
axis tight;
camlight left 
lighting none;

pause(0.2);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




