function map = element_shape_function(type,x,y,z)
%MAPFUNC  A simple transformation to map from local coordinates of unit 
%triangle and tetrahedra to global coordinates (straight edged elements!)
%MAP = mapfunc(TYPE,X,Y,Z)
%
%USE: xglobal=x1*map(1)+x2*map(2)+x3*map(3)
%     yglobal=y1*map(1)+y2*map(2)+y3*map(3)
%
%INPUT:
%1. X, Y, Z - local coordinates
%2. TYPE - string describing different element types
%   'tri3','tri6','tri10' - Coordinate mapping triangle 
%   'tet4','tet10' - Coordinate mapping tetrahedron
%
%OUTPUT
%1. MAPFUNC - vector of map functions
%
%M Crabb - 29.06.2012

if(strcmp(type(1:3),'tri'))
    map = maptri3(x,y,z);
elseif(strcmp(type(1:3),'tet'))
    map = maptet4(x,y,z);
else
    error('Incorrect number of input arguments')
end


function map = maptri3(x,y,z)
map(1) = 1-x-y;
map(2) = x;
map(3) = y;
end

function map = maptet4(x,y)
map(1) = (1-x-y)*(1-2*x-2*y);
map(2) = x*(2*x-1);
map(3) = y*(2*y-1);
map(4) = 4*x*(1-x-y);
map(5) = 4*y*(1-x-y);
map(6) = 4*x*y;
end

end
