function [vols]= check_vols(simp,vtx);
%function [vols]= check_vols(simp,vtx);
%
%Auxiliary function that calculates the volume of 
%each tetrahedron in the mesh. 
%
%
%simp = The simplices 
%vtx  = The vertices
%vols = The array of volumes.

vols = [];

for i=1:size(simp,1)
   
   this_simp = simp(i,:);
   
   t_vtx = vtx(simp(i,:)',:);
   
   t_vm = [ones(4,1),t_vtx];
   
   t_vol = abs(1/6 * det(t_vm));
   
   vols = [vols;abs(t_vol)];
   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%