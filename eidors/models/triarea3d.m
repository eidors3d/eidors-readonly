function [ta] = triarea3d(V);
%[ta] = triarea3d(V);
%
%Function that calculates the area of a triangle 
%in the 3D Cartesian space. 
%V = the 3 x 3 coordinates matrix of the points, first
%column for the xs and last for the zs.
%ta = the area of the triangle

p1 = [V(1,2) V(1,3) 1; V(2,2) V(2,3) 1; V(3,2) V(3,3) 1];      

p2 = [V(1,3) V(1,1) 1; V(2,3) V(2,1) 1; V(3,3) V(3,1) 1];

p3 = [V(1,1) V(1,2) 1; V(2,1) V(2,2) 1; V(3,1) V(3,2) 1];
 
ta = 0.5 * sqrt((det(p1))^2 + (det(p2))^2 + (det(p3))^2) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
