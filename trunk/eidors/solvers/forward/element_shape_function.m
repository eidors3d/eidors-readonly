function shape = element_shape_function(type,x,y,z)
%ELEMENT_SHAPE_FUNCTION  - Evaluation of shape functions
%values each point
%
%
%INPUT:
%1. X, Y, Z - local coordinates
%2. TYPE - string describing different element types
%   'tri3' - Linear, 2 node 1D shape function 
%   'tri 6'- Quadratic, 3 node 1D shape function
%   'tet4' - Linear, 3 node triangle shape function
%   'tet10' - Quadratic, 6 node triangle shape function
%
%OUTPUT
%1. MAPFUNC - vector of shape functions on this element
%
%M Crabb - 29.06.2012

if(strcmp(type,'tri3'))
    shape = elemshapetri3(x,y);
elseif(strcmp(type,'tri6'))
    shape = elemshapetri6(x,y);
elseif(strcmp(type,'tri10'))
    shape = elemshapetri10(x,y);
elseif(strcmp(type,'tet4'))
    shape = elemshapetet4(x,y,z);
elseif(strcmp(type,'tet10'))
    shape = elemshapetet10(x,y,z);
else
    error('Incorrect number of input arguments')
end


function map = elemshapetri3(x,y)
map(1) = 1-x-y;
map(2) = x;
map(3) = y;
end

function map = elemshapetri6(x,y)
map(1) = (1-x-y)*(1-2*x-2*y);
map(2) = x*(2*x-1);
map(3) = y*(2*y-1);
map(4) = 4*x*(1-x-y);
map(5) = 4*x*y;
map(6) = 4*y*(1-x-y);
end

function map = elemshapetri10(x,y)
map(1) = 0.5*(1-x-y)*(1-3*x-3*y)*(2-3*x-3*y);
map(2) = 0.5*x*(3*x-2)*(3*x-1);
map(3) = 0.5*y*(3*y-2)*(3*y-1);
map(4) = 4.5*x*(1-x-y)*(2-3*x-3*y);
map(5) = 4.5*x*(1-x-y)*(3*x-1);
map(6) = 4.5*x*y*(3*x-1);
map(7) = 4.5*x*y*(3*y-1);
map(8) = 4.5*y*(1-x-y)*(3*y-1);
map(9) = 4.5*y*(1-x-y)*(2-3*x-3*y);
map(10) = 27*x*y*(1-x-y);
end


function map = elemshapetet4(x,y,z)
map(1) = (1-x-y-z);
map(2) = x;
map(3) = y;
map(4) = z;
end

function map = elemshapetet10(x,y,z)
map(1) = (1-x-y-z)*(1-2*x-2*y-2*z);
map(2) = x*(2*x-1);
map(3) = y*(2*y-1);
map(4) = z*(2*z-1);
map(5) = 4*x*(1-x-y-z);
map(6) = 4*y*(1-x-y-z);
map(7) = 4*z*(1-x-y-z);
map(8) = 4*x*y;
map(9) = 4*y*z;
map(10) = 4*x*z;
end


end
