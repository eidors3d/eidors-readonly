function shape=boundary_shape_function(type,x,y)
%BOUNDSHAPEFUNC  Shape functions on the boundary in local coordiantes
%SHAPE = boundshapefunc(TYPE,X,Y)
%
%NOTES - On boundary, dimension is smaller e.g. tri6 element, boundary has 
%3 node (quadratic) and is 1D
%
%INPUT:
%1. X, Y, Z - local coordinates
%2. TYPE - string describing different element typeswhat kind of element
%   'tri3' - Linear, 2 node 1D shape function 
%   'tri 6'- Quadratic, 3 node 1D shape function
%   'tet4' - Linear, 3 node triangle shape function
%   'tet10' - Quadratic, 6 node triangle shape function
%
%OUTPUT
%1. SHAPE - vector of shape functions on this segment of boundary
%
%M Crabb - 29.06.2012

if(strcmp(type,'tri3'))
    shape = boundshapetri3(x,y);
elseif(strcmp(type,'tri6'))
    shape = boundshapetri6(x,y);
elseif(strcmp(type,'tri10'))
    shape = boundshapetri10(x,y);
elseif(strcmp(type,'tet4'))
    shape = boundshapetet4(x,y);
elseif(strcmp(type,'tet10'))
    shape = boundshapetet10(x,y);
else
    error('Incorrect number of input arguments')
end



%1D Shape functions : Equidistant placed in region [-1,1]

function shape = boundshapetri3(x,y)
shape(1) = 0.5*(1-x);
shape(2) = 0.5*(1+x);
end

%Ordered correctly (vertex basis functions are first : (1,3,2) order)
function shape = boundshapetri6(x,y)
shape(1) = 0.5*x*(x-1);
shape(2) = 0.5*x*(x+1);
shape(3) = (1-x)*(1+x);
end

%Ordered correctly (vertex basis functions are first : (1,3,4,2) order)
function shape = boundshapetri10(x,y)
shape(1) = 0.0625*(1-x)*(3*x-1)*(3*x+1);
shape(2) = 0.0625*(1+x)*(3*x-1)*(3*x+1);
shape(3) = -0.5625*(1+x)*(3*x-1)*(1-x);
shape(4) = 0.5625*(1+x)*(3*x+1)*(1-x);
end





%2D Shape functions : Place on the unit triangle [0,1]*[0,1] axes

function shape = boundshapetet4(x,y)
shape(1) = 1-x-y;
shape(2) = x;
shape(3) = y;
end

%Ordered correctly (i.e. vertex basis function (first 3) are first)
function shape = boundshapetet10(x,y)
shape(1) = (1-x-y)*(1-2*x-2*y);
shape(2) = x*(2*x-1);
shape(3) = y*(2*y-1);
shape(4) = 4*x*(1-x-y);
shape(5) = 4*y*(1-x-y);
shape(6) = 4*x*y;
end

end
