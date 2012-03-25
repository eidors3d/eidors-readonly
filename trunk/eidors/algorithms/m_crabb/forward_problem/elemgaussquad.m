 function [w,x,y,z] = elemgaussquad(type)
%ELEMGAUSSQUAD  
%Weights and coordinate for Gauss integration over reference element as a
%function of the element type. Specifically, we have rules that integrate
%products of gradients of shape functions.
%[w,x,y,z] = elemgaussquad(TYPE)
%
%IMPORTANT NOTE: Integrals of f on reference element
%       I = int(f(x,y,z)) dV =  w_{i}*f(x_{i},y_{i},z_{i}
%Reference triangle defined as x,y,z - ([0,1],[0,1],[0,1])
%
%INPUT:
%1. TYPE - string describing different element types
%   'tri3' - Linear, 3 node triangle 
%   'tri6'- Quadratic, 6 node triangle
%   'tri10' - Cubic, 10 node triangle
%   'tet4' - Linear, 4 node tetrahedral
%   'tet10' - Quadratic, 10 node tetrahedral
%
%OUTPUT
%1. w - the weights matrix of shape function derivatives size(ndim,nshape)
%2. x,y,z - evaluation points of the function f

if(strcmp(type,'tri3'))
    [w,x,y,z] = elemgaussquadtri3;
elseif(strcmp(type,'tri6'))
    [w,x,y,z] = elemgaussquadtri6;
elseif(strcmp(type,'tri10'))
    [w,x,y,z] = elemgaussquadtri10;
elseif(strcmp(type,'tet4'))
    [w,x,y,z] = elemgaussquadtet4;
elseif(strcmp(type,'tet10'))
    [w,x,y,z] = elemgaussquadtet10;
elseif(strcmp(type,'order5'))
    [w,x,y,z]=elemgaussquadtripolyorder5;
elseif(strcmp(type,'order7'))
    [w,x,y,z]=elemgaussquadtripolyorder7;
else
    error('Integration on this element type not recognised');
end



%TRIANGLE - USING [0,1]*[0,1] UNIT REFERENCE TRIANGLE

%Order 0 polynomial (f=a)
function [w,x,y,z] = elemgaussquadtri3
    %Integration scheme integrates linear exactly on unit triangle
    w(1)=1.0/2.0; x(1)=1.0/3.0; y(1)=1.0/3.0; z(1)=0.0;
end

%Order 2 polynomial (f=a+b*x+c*y+d*x*x+e*x*y+f*y*y)
function [w,x,y,z] = elemgaussquadtri6
    %Integration scheme integrate quadratic exactly on unit triangle
    w(1)=1.0/6.0; x(1)=2.0/3.0; y(1)=1.0/6.0; z(1)=0.0;
    w(2)=1.0/6.0; x(2)=1.0/6.0; y(2)=1.0/6.0; z(2)=0.0;
    w(3)=1.0/6.0; x(3)=1.0/6.0; y(3)=2.0/3.0; z(3)=0.0;   
end

%Order 4 polynomial (f=a+b*x+c*y+d*x*x+e*x*y+f*y*y+g*x*x*x+h*x*x*y+i*x*y*y+j*y*y*y)
function [w,x,y,z] = elemgaussquadtri10
%Integration scheme integrates quartic exactly on unit triangle
    w(1)=0.5*0.2233815896780; x(1)=0.1081030181681; y(1)=0.4459484909160; z(1)=0.0;
    w(2)=0.5*0.2233815896780; x(2)=0.4459484909160; y(2)=0.1081030181681; z(2)=0.0;
    w(3)=0.5*0.2233815896780; x(3)=0.4459484909160; y(3)=0.4459484909160; z(3)=0.0;
    w(4)=0.5*0.1099517436553; x(4)=0.8168475729805; y(4)=0.0915762135098; z(4)=0.0;
    w(5)=0.5*0.1099517436553; x(5)=0.0915762135098; y(5)=0.0915762135098; z(5)=0.0;
    w(6)=0.5*0.1099517436553; x(6)=0.0915762135098; y(6)=0.8168475729805; z(6)=0.0;
end

function [w,x,y,z]=elemgaussquadtripolyorder5
    %Integration scheme integrates quintic exactly on unit triangle
    w(1)=0.5*0.1259391805448; x(1)=0.1012865073235; y(1)=0.1012865073235; z(1)=0.0;
    w(2)=0.5*0.1259391805448; x(2)=0.7974269853531; y(2)=0.1012865073235; z(2)=0.0;
    w(3)=0.5*0.1259391805448; x(3)=0.1012865073235; y(3)=0.7974269853531; z(3)=0.0;
    w(4)=0.5*0.1323941527885; x(4)=0.4701420641051; y(4)=0.0597158717898; z(4)=0.0;
    w(5)=0.5*0.1323941527885; x(5)=0.4701420641051; y(5)=0.4701420641051; z(5)=0.0;
    w(6)=0.5*0.1323941527885; x(6)=0.0597158717898; y(6)=0.4701420641051; z(6)=0.0;
    w(7)=0.5*0.225;           x(7)=0.3333333333333; y(7)=0.3333333333333; z(7)=0.0;
end

function [w,x,y,z]=elemgaussquadtripolyorder7
    %Integration scheme integrate heptic exactly on unit triangle
    w(1)=0.5*0.0533472356088;x(1)=0.0651301029022; y(1)=0.0651301029022; z(1)=0.0;
    w(2)=0.5*0.0533472356088;x(2)=0.8697397941956; y(2)=0.0651301029022; z(2)=0.0;
    w(3)=0.5*0.0533472356088;x(3)=0.0651301029022; y(3)=0.8697397941956; z(3)=0.0;
    w(4)=0.5*0.0771137608903;x(4)=0.3128654960049; y(4)=0.0486903154253; z(4)=0.0;
    w(5)=0.5*0.0771137608903;x(5)=0.6384441885698; y(5)=0.3128654960049; z(5)=0.0;
    w(6)=0.5*0.0771137608903;x(6)=0.0486903154253; y(6)=0.6384441885698; z(6)=0.0;
    w(7)=0.5*0.0771137608903;x(7)=0.6384441885698; y(7)=0.0486903154253; z(7)=0.0;
    w(8)=0.5*0.0771137608903;x(8)=0.3128654960049; y(8)=0.6384441885698; z(8)=0.0;
    w(9)=0.5*0.0771137608903;x(9)=0.0486903154253; y(9)=0.3128654960049; z(9)=0.0;
    w(10)=0.5*0.1756152574332;x(10)=0.2603459660790; y(10)=0.2603459660790; z(10)=0.0;
    w(11)=0.5*0.1756152574332;x(11)=0.4793080678419; y(11)=0.2603459660790; z(11)=0.0;
    w(12)=0.5*0.1756152574332;x(12)=0.2603459660790; y(12)=0.4793080678419; z(12)=0.0;
    w(13)=-0.5*0.1495700444677;x(13)=0.3333333333333; y(13)=0.3333333333333; z(13)=0.0;
end



%TETRAHEDRON - USING [0,1]*[0,1]*[0,1] UNIT REFERENCE TETRAHEDRON

%Order 0 polynomial (f=a)
function [w,x,y,z] = elemgaussquadtet4
    %Integration scheme integrates exactly linear on unit tetrahedron
    w(1)=1.0/6.0; x(1)=1.0/4.0; y(1)=1.0/4.0; z(1)=1.0/4.0;
end

%Order 2 polynomial (f=a+b*x+c*y+d*z+e*x*x+f*y*y+g*z*z+h*x*y+i*y*z+j*x*z)
function [w,x,y,z] = elemgaussquadtet10
    %Integration scheme integrates exactly quadratic on unit tetrahedron
    a=(5+3*sqrt(5))/20; b=(5-sqrt(5))/20;
    w(1)=1.0/24.0; x(1)=a; y(1)=b; z(1)=b;
    w(2)=1.0/24.0; x(2)=b; y(2)=a; z(2)=b;
    w(3)=1.0/24.0; x(3)=b; y(3)=b; z(3)=a;
    w(4)=1.0/24.0; x(4)=b; y(4)=b; z(4)=b;
end

end
