function [w,x,y]=boundary_gauss_points(type)
%BOUNDGAUSSQUAD  Weights and coordinate for Gauss integration over 
%over triangles as a function of the element type
%[w,x,y,z] = boundgaussquad(TYPE,X,Y,Z)
%
%IMPORTANT NOTE: Integrals of f on reference element
%       I = int(f(x,y,z)) dV =  w_{i}*f(x_{i},y_{i},z_{i}
%1D boundary - Reference line defined as x - [-1,1]
%2D boundary - Reference triangle defined as x,y - ([0,1],[0,1])
%
%On boundary, dimension is smaller e.g. tri6 element, boundary has 
%node (quadratic) and is 1D
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
%1. w - the weights matrix of shape function derivatives size(ndim,nshape)
%2. x,y,z - evaluation points of the function f
%
%M Crabb - 29.06.2012

%Test what type of element we have
if(strcmp(type,'tri3'))
    [w,x,y] = boundgaussquadtri3;
elseif(strcmp(type,'tri6'))
    [w,x,y] = boundgaussquadtri6;
elseif(strcmp(type,'tri10'))
    [w,x,y] = boundgaussquadtri10;
elseif(strcmp(type,'order9'))
    [w,x,y]=boundgaussquadorder9;
elseif(strcmp(type,'tet4'))
    [w,x,y] = boundgaussquadtet4;
elseif(strcmp(type,'tet10'))
    [w,x,y] = boundgaussquadtet10;
else
    error('Integration on this element type not recognised');
end

%LINEAR - USING [-1,1]
%Order 2 polynomial on triangle edge (need 2 point Gauss rule )
function [w,x,y] = boundgaussquadtri3
    %Integration scheme integrates quadratic exactly on [-1,1]
    w(1)=1.0; x(1)=1.0/sqrt(3.0);     y(1)=0.0;
    w(2)=1.0; x(2)=-1.0/sqrt(3.0);    y(2)=0.0;
end

%Order 4 polynomial on triangle edge (need 3 point Gauss ruler)
function [w,x,y] = boundgaussquadtri6
    %Integration scheme integrates quartic exactly on [-1,1]
    w1=sqrt(3.0/5.0);
    w(1)=8.0/9.0; x(1)=0.0; y(1)=0.0;
    w(2)=5.0/9.0; x(2)=w1;  y(2)=0.0;
    w(3)=5.0/9.0; x(3)=-w1; y(3)=0.0;
end

%Order 6 polynomial on triangle edge (need 4 point Gauss rule)
function [w,x,y] = boundgaussquadtri10
    %Integration scheme integrates hexic exactly on [-1,1]
    w1=(18.0+sqrt(30.0))/36.0; w2=(18.0-sqrt(30.0))/36.0;
    p1=sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0);
    p2=sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0);
    w(1)=w1;   x(1)=p1;    y(1)=0.0;
    w(2)=w1;   x(2)=-p1;   y(2)=0.0;
    w(3)=w2;   x(3)=p2;    y(3)=0.0;
    w(4)=w2;   x(4)=-p2;   y(4)=0.0;
end

%5 point rule accurate for polynomials up to degree 9
function [w,x,y] = boundgaussquadorder9
    %Integration scheme integrates hexic exactly on [-1,1]
    w1=128/225; 
    w2=(322.0 + 13.0*sqrt(70.0))/900.0;
    w3=(322.0 - 13.0*sqrt(70.0))/900.0;
    p1=sqrt( 5.0 - 2.0*sqrt(10.0/7.0) )/3.0;
    p2=sqrt( 5.0 + 2.0*sqrt(10.0/7.0) )/3.0;
    w(1)=w1;   x(1)=0;    y(1)=0.0;
    w(2)=w2;   x(2)=p1;   y(2)=0.0;
    w(3)=w2;   x(3)=-p1;    y(3)=0.0;
    w(4)=w3;   x(4)=p2;   y(4)=0.0;
    w(5)=w3;   x(5)=-p2;   y(5)=0.0;
end




%TRIANGLE - USING [0,1]*[0,1]
%Order 2 polynomial on tetrahedron face (quadratic triangle rule)
function [w,x,y] = boundgaussquadtet4
    %Integration scheme integrates quadratic exactly on unit triangle
    w(1)=1.0/6.0; x(1)=2.0/3.0; y(1)=1.0/6.0;
    w(2)=1.0/6.0; x(2)=1.0/6.0; y(2)=1.0/6.0;
    w(3)=1.0/6.0; x(3)=1.0/6.0; y(3)=2.0/3.0;
end

%Order 4 polynomial on tetrahedron face (quartic triangle rule) 
function [w,x,y] = boundgaussquadtet10
    %Integration scheme integrates quartic exactly on unit triangle
    w(1)=0.5*0.2233815896780; x(1)=0.1081030181681; y(1)=0.4459484909160;
    w(2)=0.5*0.2233815896780; x(2)=0.4459484909160; y(2)=0.1081030181681;
    w(3)=0.5*0.2233815896780; x(3)=0.4459484909160; y(3)=0.4459484909160;
    w(4)=0.5*0.1099517436553; x(4)=0.8168475729805; y(4)=0.0915762135098;
    w(5)=0.5*0.1099517436553; x(5)=0.0915762135098; y(5)=0.0915762135098;
    w(6)=0.5*0.1099517436553; x(6)=0.0915762135098; y(6)=0.8168475729805;    
end

end
