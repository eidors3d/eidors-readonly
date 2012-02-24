function shape = elemshapefunc(type,x,y,z)
%ELEMSHAPEFUNC  Shape functions of elements in local coordiantes
%SHAPE = elemshapefunc(TYPE,X,Y,Z)
%
%INPUT:
%1. X, Y, Z - local coordinates
%2. TYPE - string describing different element typeswhat kind of element
%   'tri3' - Linear, 3 node triangle 
%   'tri 6'- Quadratic, 6 node triangle
%   'tet4' - Linear, 4 node tetrahedral
%   'tet10' - Quadratic, 10 node tetrahedral
%
%OUTPUT
%1. SHAPE - vector of shape functions

if(strcmp(type,'tri3'))
    shape = shapetri3(x,y,z);
elseif(strcmp(type,'tri6'))
    shape = shapetri6(x,y,z);
elseif(strcmp(type,'tri10'))
    shape = shapetri10(x,y,z);    
elseif(strcmp(type,'tet4'))
    shape = shapetet4(x,y,z);
elseif(strcmp(type,'tet10'))
    shape = shapetet10(x,y,z);
else
    error('Incorrect number of input arguments')
end


function shape = shapetri3(x,y,z)
shape(1) = 1-x-y;
shape(2) = x;
shape(3) = y;
end

function shape = shapetri6(x,y,z)
shape(1) = (1-x-y)*(1-2*x-2*y);
shape(2) = x*(2*x-1);
shape(3) = y*(2*y-1);
shape(4) = 4*x*(1-x-y);
shape(5) = 4*x*y;
shape(6) = 4*y*(1-x-y);
end

function shape = shapetri10(x,y,z)
shape(1) = 0.5*(1-x-y)*(1-3*x-3*y)*(2-3*x-3*y);
shape(2) = 0.5*x*(3*x-2)*(3*x-1); 
shape(3) = 0.5*y*(3*y-2)*(3*y-1);
shape(4) = 4.5*x*(1-x-y)*(2-3*x-3*y);
shape(5) = 4.5*x*(1-x-y)*(3*x-1);
shape(6) = 4.5*x*y*(3*x-1);
shape(7) = 4.5*x*y*(3*y-1);
shape(8) = 4.5*y*(1-x-y)*(3*y-1);
shape(9) = 4.5*y*(1-x-y)*(2-3*x-3*y);
shape(10) = 27*x*y*(1-x-y);
end

function shape = shapetet4(x,y,z)
shape(1) = 1-x-y-z;
shape(2) = x;
shape(3) = y;
shape(4) = z;
end

function shape = shapetet10(x,y,z)
shape(1) = (1-x-y-z)*(1-2*x-2*y-2*z);
shape(2) = x*(2*x-1);
shape(3) = y*(2*y-1);
shape(4) = z*(2*z-1);
shape(5) = 4*x*(1-x-y-z);
shape(6) = 4*y*(1-x-y-z);
shape(7) = 4*z*(1-x-y-z);
shape(8) = 4*x*y;
shape(9) = 4*y*z;
shape(10) = 4*z*x;
end

end
