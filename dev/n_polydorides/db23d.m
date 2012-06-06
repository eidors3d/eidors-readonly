function [dd] = db23d(x1,y1,z1,x2,y2,z2);
%function [dd] = db23d(x1,y1,z1,x2,y2,z2);
%
%Auxiliary function that caclulates the distance between 
%two points or two sets of points in 3D
%
%
%
%(x1,y1,z1) = The coordinates of the first point(s) in 3D
%(x2,y2,z2) = The coordinates of the second point(s) in 3D
%dd         = Their distance(s)

dd = sqrt((x2 - x1).^2 + (y2 - y1).^2 + (z2 - z1).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%