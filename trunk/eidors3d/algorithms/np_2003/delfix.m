function [vtx_n,simp_n] = delfix(vtx,simp)
%function [vtx_n,simp_n] = delfix(vtx,simp)
%
% Auxiliary function to remove the zero area faces
% produced by Matlab's delaunay triangulation
%
%
%
%vtx  = The vertices matrix
%simp = The simplices matrix


simp_n = [];
tri_a = [];

for kk=1:length(simp)
   
   this_tri = simp(kk,:);
   
   xa = vtx(this_tri(1),1); ya = vtx(this_tri(1),2);
   xb = vtx(this_tri(2),1); yb = vtx(this_tri(2),2);
   xc = vtx(this_tri(3),1); yc = vtx(this_tri(3),2);
   
   tria = polyarea([xa;xb;xc],[ya;yb;yc]);
   tri_a = [tri_a ; tria];
   
   if tria > 0.00000000001       
      simp_n = [simp_n;this_tri];
   end
   
end

vtx_n = vtx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
