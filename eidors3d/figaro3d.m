function [fc] = figaro3d(srf,vtx,simp,fc,BB,h);
%function [fc] = figaro3d(srf,vtx,simp,fc,BB,h);
%
%This function plots the solution as a 3D object crossed at the plane z=h 
%within the 3D boundaries of the volume
%
%
%
%h    = is the height of the desired slice, max(vtx(:,3))>= h >= min(vtx(:,3))
%srf  = The boundary surfaces (faces)
%simp = The simplices matrix
%fc   = Plotting utility matrix, calculated through slicer_plot. 
%       This is calculated in slicer_plot.m 
%BB   = The calculated inverse solution

srf_m=[];

for i=1:size(srf,1)
   
   if vtx(srf(i,1),3) < h & vtx(srf(i,2),3) < h & vtx(srf(i,3),3) < h 
      %Adjust to <= if a coarse mesh is used.
      
      srf_m = [srf_m;srf(i,:)];
   end
end

figure;

p = trisurf(srf_m,vtx(:,1),vtx(:,2),vtx(:,3),ones(size(vtx,1),1));
set(p,'FaceColor',[0.6 1 0.5],'EdgeColor','none');
camlight left 
lighting flat;
hold on;
[fc] = slicer_plot(h,BB,vtx,simp,fc);
%or change to
%[fc] = slicer_plot_n(h,BB,vtx,simp,fc);
daspect([1 1 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
