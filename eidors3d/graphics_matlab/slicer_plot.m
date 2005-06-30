function [fc] = slicer_plot(h,BB,vtx,simp,fc);
%function [fc] = slicer_plot(h,BB,vtx,simp,fc);
%
%This function plots a 2D slice of the 3D solution vector BB at z=h, by
%extracting a smoother, interpolated, node-wise solution from the caclulated
%elent-wise inverse solution.
%
%
%
%h    = The height of the desired solution, max(vtx(:,3))>= h >= min(vtx(:,3)).
%BB   = The caclulated inverse solution
%vtx  = The vertices matrix
%simp = The simplices matrix
%fc   = The edges of the mesh. This is a 2 column matrix required for the 3D plotting. 
%       fc may take some time to be calculated so it is a good idea to save it and use 
%       it thereafter. Initially use [fc] = slicer_plot(h,BB,vtx,simp); to plot the slide 
%       and calculate fc.


if nargin < 5 || (nargin == 6 && isempty(fc))
fc = [];
%Develop the faces           

for f=1:size(simp,1)
   
   fc1 = sort([simp(f,1),simp(f,2)]);
   fc2 = sort([simp(f,1),simp(f,3)]);
   fc3 = sort([simp(f,1),simp(f,4)]);
   fc4 = sort([simp(f,2),simp(f,3)]);
   fc5 = sort([simp(f,2),simp(f,4)]);
   fc6 = sort([simp(f,3),simp(f,4)]);
   
   fc = [fc;fc1;fc2;fc3;fc4;fc5;fc6];
   
end

fc = unique(fc,'rows');

end


%Extract the nodal solution CC
[KK] = solution_ext(BB,vtx,simp);

%Need to rescale that to match BB.

ranges = (max(BB)-min(BB))/(max(KK)-min(KK));

CC = zeros(size(KK));

for yu=1:length(KK)
   CC(yu) = ranges * KK(yu);
end


vtxp = []; %Nodes created for the plane
vap = []; %Value of the node in vtxp

for j=1:size(fc,1)
   
   this_ph = fc(j,:); %[nodeA nodeB]
   
   if max(vtx(this_ph(1),3),vtx(this_ph(2),3))> h & ...
         min(vtx(this_ph(1),3),vtx(this_ph(2),3))<= h 		
     
  %If the face is crossed through by the plane then 
  %create a plotable node on the plane.
  
  Pa = this_ph(1); Pb = this_ph(2);
  
  xa = vtx(Pa,1); ya = vtx(Pa,2); za = vtx(Pa,3);
  xb = vtx(Pb,1); yb = vtx(Pb,2); zb = vtx(Pb,3);
  
  
  xn = xa + (h-za)*(xb-xa)/(zb-za);
  yn = ya + (h-za)*(yb-ya)/(zb-za);
  
  vtxp = [vtxp;[xn,yn]];
  
  %Specifing the value of the new node
  
  [la] = db23d(xa,ya,za,xn,yn,h);
  
  [lb] = db23d(xb,yb,zb,xn,yn,h);
  
  vv = (CC(Pa)*(lb)/(la + lb)) + (CC(Pb)*(la)/(la + lb));
  
  vap = [vap;vv];
  
end %if

end %for

tri = delaunay(vtxp(:,1),vtxp(:,2));

[vtxp,tri] = delfix(vtxp,tri);

X = vtxp(:,1);
Y = vtxp(:,2);
Z = h*ones(length(vtxp),1);
trisurf(tri,X,Y,Z,vap);
set(gca,'LineWidth',1e-6);
shading interp;  %%%%%%If the solution is too smooth and you can't see it comment this line.
axis([min(vtx(:,1)), max(vtx(:,1)), ...
      min(vtx(:,2)), max(vtx(:,2)), ...
      min(vtx(:,3)), max(vtx(:,3))]);

%You may want to change this setting !!  
caxis([min(BB), max(BB)]);

axis square;
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  



