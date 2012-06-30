function dshape = element_d_shape_function(type,x,y,z)
%DELEMSHAPEFUNC  
%Derivative of shape functions of elements in local coordiantes. The basis 
%functions are for an anti-clockwise vertex arrangement triangle (or in 3D
%anti-clockwise to the point of looking INTO triangle from a vertex)
%DSHAPE = dshapefunc(TYPE,X,Y,Z)
%
%INPUT:
%1. X, Y, Z - local coordinates
%2. TYPE - string describing different element typeswhat kind of element
%   'tri3'  - Linear, 3 node triangle 
%   'tri6'  - Quadratic, 6 node triangle
%   'tri10' - Cubic, 10 node triangle
%   'tet4'  - Linear, 4 node tetrahedral
%   'tet10' - Quadratic, 10 node tetrahedral
%
%OUTPUT
%1. DELEMSHAPEFUNC - matrix of shape function derivatives size(ndim,nshape)
%              - using convention : dshape(i,j) = d shape_{j} / d x_{i}
%
%M Crabb - 29.06.2012

if(strcmp(type,'tri3'))
    dshape = delemshapetri3(x,y,z);
elseif(strcmp(type,'tri6'))
    dshape = delemshapetri6(x,y,z);
elseif(strcmp(type,'tri10'))
    dshape = delemshapetri10(x,y,z);
elseif(strcmp(type,'tet4'))
    dshape = delemshapetet4(x,y,z);
elseif(strcmp(type,'tet10'))
    dshape = delemshapetet10(x,y,z);
else
    error('Incorrect number of input arguments')
end

%TRIANGLE - USING [0,1]*[0,1] UNIT REFERENCE TRIANGLE

function dshape = delemshapetri3(x,y,z)
dshape(1,1) = -1; dshape(1,2) = 1; dshape(1,3) = 0;
dshape(2,1) = -1; dshape(2,2) = 0; dshape(2,3) = 1;
end

function dshape = delemshapetri6(x,y,z)
dshape(1,1) = -3 + 4*x + 4*y; dshape(2,1) = -3 + 4*x + 4*y;
dshape(1,2) = 4*x - 1;        dshape(2,2) = 0;
dshape(1,3) = 0;              dshape(2,3) = 4*y-1;
dshape(1,4) = 4 - 8*x - 4*y;  dshape(2,4) = -4*x;
dshape(1,5) = 4*y;            dshape(2,5) = 4*x; 
dshape(1,6) = -4*y;           dshape(2,6) = 4 - 4*x - 8*y;
end

function dshape=delemshapetri10(x,y,z)
dshape(1,1) = 0.5*(36*x+36*y-27*x*x-27*y*y-54*x*y-11); dshape(2,1) = 0.5*(36*x+36*y-27*x*x-27*y*y-54*x*y-11);
dshape(1,2) = 13.5*x*x-9*x+1;                          dshape(2,2) = 0;
dshape(1,3) = 0;                                       dshape(2,3) = 13.5*y*y-9*y+1;
dshape(1,4) = 4.5*(2-10*x+9*x*x+12*x*y-5*y+3*y*y);     dshape(2,4) = 4.5*(6*x*x-5*x+6*x*y);
dshape(1,5) = 4.5*(y-6*x*y-9*x*x+8*x-1);               dshape(2,5) = 4.5*x*(1-3*x);
dshape(1,6) = 27*x*y-4.5*y;                            dshape(2,6) = 4.5*x*(3*x-1);
dshape(1,7) = 4.5*y*(3*y-1);                           dshape(2,7) = 27*x*y-4.5*x;
dshape(1,8) = 4.5*y*(1-3*y);                           dshape(2,8) = 4.5*(x-6*x*y-9*y*y+8*y-1);
dshape(1,9) = 4.5*(6*x*y-5*y+6*y*y);                   dshape(2,9) = 4.5*(2-5*x+3*x*x+12*x*y-10*y+9*y*y);
dshape(1,10) = 27*y-54*x*y-27*y*y;                     dshape(2,10) = 27*x-54*x*y-27*x*x;
end






%TETRAHEDRON - USING [0,1]*[0,1]*[0,1] UNIT REFERENCE TETRAHEDRON

function dshape = delemshapetet4(x,y,z)
dshape(1,1) = -1; dshape(2,1) = -1; dshape(3,1) = -1;
dshape(1,2) = 1;  dshape(2,2) = 0;  dshape(3,2) = 0;
dshape(1,3) = 0;  dshape(2,3) = 1;  dshape(3,3) = 0;
dshape(1,4) = 0;  dshape(2,4) = 0;  dshape(3,4) = 1;
end

function dshape = delemshapetet10(x,y,z)
dshape(1,1) = -3 + 4*x + 4*y + 4*z; dshape(2,1) = -3 + 4*x + 4*y + 4*z; dshape(3,1) = -3 + 4*x + 4*y + 4*z;
dshape(1,2) = 4*x - 1;              dshape(2,2) = 0;                    dshape(3,2) = 0; 
dshape(1,3) = 0;                    dshape(2,3) = 4*y-1;                dshape(3,3) = 0;
dshape(1,4) = 0;                    dshape(2,4) = 0;                    dshape(3,4) = 4*z - 1;
dshape(1,5) = 4 - 8*x - 4*y - 4*z;  dshape(2,5) = -4*x;                 dshape(3,5) = -4*x;
dshape(1,6) = -4*y;                 dshape(2,6) = 4 - 4*x - 8*y - 4*z;  dshape(3,6) = -4*y;
dshape(1,7)= -4*z;                  dshape(2,7)= -4*z;                  dshape(3,7) = 4 - 4*x - 4*y - 8*z;
dshape(1,8) = 4*y;                  dshape(2,8) = 4*x;                  dshape(3,8) = 0;
dshape(1,9) = 0;                    dshape(2,9) = 4*z;                  dshape(3,9) = 4*y;
dshape(1,10) = 4*z;                 dshape(2,10) = 0;                   dshape(3,10) = 4*x;
end

end
