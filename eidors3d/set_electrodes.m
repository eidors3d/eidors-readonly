function [elec_face,sels,cnts,VV] = set_electrodes(vtx,srf,elec_face,sels,cnts,VV);
%function [elec_face,sels,cnts,VV] = set_electrodes(vtx,srf,elec_face,sels,cnts,VV);
%
%This function must be called recursively to sellect boundary faces building up 
%the elec_face. You will need to reshape this matrix appropriately to get the elec 
%matrix, depending on how many faces there are in each electrode.
%
%
%
%vtx      = The vertices matrix
%srf      = The boundary surfaces
%elec_face = A 3 column matrix holding the boundary faces to be used for constructing 
%           the electrodes
%sels     = The indices in srf matrix of the sellected surfaces
%cnts     = The coordinates of the center of each triangular boundary surface.
%VV       = The last viewing angle.
%
%Call this function as follows
%
%[elec_face,sels,cnts,VV] = set_electrodes(vtx,srf,[],[],[],[30 60]);
%ONLY FOR THE FIRST TIME, for the first row of elec_face.
%and then ...
%[elec_face,sels,cnts,VV] = set_electrodes(vtx,srf,elec_face,sels,cnts,VV);
%Once you have gathered the faces to be assigned as electrodes,  
%provided that you have (the same) N number of faces per electrode and these are indexed in 
%sequence inside "elec_face", then use something like
%"elec = reshape(elec_face',(N*3),number of electrodes)';"
%In this case note that "elec_face = srf(sels,:);" should be valid.


if size(cnts,1) == 0


cnts = []; % Vector containing the geometric center of each triangular 
			  % surface in x,y,z coordinates.
 
for i=1:size(srf,1)
   
   a = srf(i,1);
   b = srf(i,2);
   c = srf(i,3);
   
   ccnx = (vtx(a,1) + vtx(b,1) + vtx(c,1))/3;
   ccny = (vtx(a,2) + vtx(b,2) + vtx(c,2))/3;
   ccnz = (vtx(a,3) + vtx(b,3) + vtx(c,3))/3;
   
   ccn = [ccnx,ccny,ccnz];
   
   cnts = [cnts; ccn];
end
end

%Plot the surface

trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
grid off
hidden off %%%% If you have viewing problems comment this line
view(VV);


%and any previous patches
if ~isempty(sels) == 1 
for y=1:size(sels)
      paint_electrodes(sels(y),srf,vtx);
end
end

disp('Rotate and click on the figure to locate electrode and then press ENTER');
disp('Try to aim near the centre of the triangular face from a small angle');

pause; 

VV = get(gca,'View');

[sel] = laserbeam(vtx,srf,cnts);

paint_electrodes(sel,srf,vtx);
   

%Now confirmation about sel is required.

button = questdlg('Was the positioning OK?','Electrode confirmation','Yes','No','Help','No');
if strcmp(button,'Yes')
   disp('Creating electrode');
elseif strcmp(button,'No')
   disp('Canceling electrode');
   elseif strcmp(button,'Help')
   disp('Sorry, no help available')
end


if strcmp(button,'Yes')
elec_face = [elec_face; srf(sel,:)];
sels = [sels;sel];
end

   
if strcmp(button,'No')  
   
elec_face = elec_face;
   
%Plot the surface
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
grid off
hidden off %%%% If you have viewing problems, comment this line

   for u=1:size(sels)
       paint_electrodes(sels(u),srf,vtx);
   end
  
end

view(VV);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%