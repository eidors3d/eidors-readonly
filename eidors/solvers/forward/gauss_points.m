function [w,x,y,z]=gauss_points(dim,order)
%GAUSS_POINTS
%Weights and coordinate for Gauss integration over reference element as a
%function of the element type. Specifically, we have rules that integrate
%products of gradients of shape functions.
%[w,x,y,z] = gauss_points(DIM,ORDER)
%
%NOTE: Integrals of f on specific reference element
%       I = int(f(x,y,z)) dV =  w_{i}*f(x_{i},y_{i},z_{i}
%Reference triangle dependent on dimension
%   DIM=1 - Reference = [-1,1]
%   DIM=2 - Reference = [0,1]x[0,1]     
%   DIM=3 - Reference = [0,1]x[0,1]x[0,1]
%
%INPUT:
%1. DIM - Dimension of integration rule - 1,2,3
%2. ORDER - Order integration exact for - rules below are written up
%   DIM = 1, ORDER=2,4,6,9
%   DIM = 2, ORDER=0,2,4,5,7
%   DIM = 3, ORDER=0,2,4
%OUTPUT
%1. w - the weights matrix of shape function derivatives size(ndim,nshape)
%2. x,y,z - evaluation points of the function f
%
%M Crabb - 29.06.2012


%Switch between cases of dimension and polynomial order
if(dim==1) 
    if(order==2)
        [w,x,y,z]=dim_1_order_2;
    elseif(order==4)
        [w,x,y,z]=dim_1_order_4;
    elseif(order==6)
        [w,x,y,z]=dim_1_order_6;        
    elseif(order==9)
        [w,x,y,z]=dim_1_order_9;
    else
       error('Integration order not recognised for dimension'); 
    end
elseif(dim==2)
    if(order==0)
        [w,x,y,z]=dim_2_order_0;       
    elseif(order==2)
        [w,x,y,z]=dim_2_order_2;        
    elseif(order==4)
        [w,x,y,z]=dim_2_order_4;        
    elseif(order==5)
        [w,x,y,z]=dim_2_order_5;        
    elseif(order==7)
        [w,x,y,z]=dim_2_order_7;
    else
        error('Integration order not recognised for dimension');
    end            
elseif(dim==3)
    if(order==0)
        [w,x,y,z]=dim_3_order_0;        
    elseif(order==2)
        [w,x,y,z]=dim_3_order_2;        
    elseif(order==4)
        [w,x,y,z]=dim_3_order_4;    
    else
        error('Integration order not recognised for dimension');
    end
    
else
    error('Dimension not recognised');
end

%1 DIMENSIONAL INTEGRATION RULES
%Order 2 polynomial on triangle edge (need 2 point Gauss rule )
function [w,x,y,z] = dim_1_order_2
    %Integration scheme integrates quadratic exactly on [-1,1]
    w(1)=1.0; x(1)=1.0/sqrt(3.0);     y(1)=0.0; z(1)=0.0;
    w(2)=1.0; x(2)=-1.0/sqrt(3.0);    y(2)=0.0; z(2)=0.0;
end

%Order 4 polynomial on triangle edge (need 3 point Gauss ruler)
function [w,x,y,z] = dim_1_order_4
    %Integration scheme integrates quartic exactly on [-1,1]
    w1=sqrt(3.0/5.0);
    w(1)=8.0/9.0; x(1)=0.0; y(1)=0.0; z(1)=0.0;
    w(2)=5.0/9.0; x(2)=w1;  y(2)=0.0; z(2)=0.0;
    w(3)=5.0/9.0; x(3)=-w1; y(3)=0.0; z(3)=0.0;
end

%Order 6 polynomial on triangle edge (need 4 point Gauss rule)
function [w,x,y,z] = dim_1_order_6
    %Integration scheme integrates hexic exactly on [-1,1]
    w1=(18.0+sqrt(30.0))/36.0; w2=(18.0-sqrt(30.0))/36.0;
    p1=sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0);
    p2=sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0);
    w(1)=w1;   x(1)=p1;    y(1)=0.0; z(1)=0.0;
    w(2)=w1;   x(2)=-p1;   y(2)=0.0; z(2)=0.0;
    w(3)=w2;   x(3)=p2;    y(3)=0.0; z(3)=0.0;
    w(4)=w2;   x(4)=-p2;   y(4)=0.0; z(4)=0.0;
end

%Order 9 polynomial on triangle edge (need 5 point Gauss rule)
function [w,x,y,z] = dim_1_order_9
    %Integration scheme integrates nonic exactly on [-1,1]
    w1=128/225; 
    w2=(322.0 + 13.0*sqrt(70.0))/900.0;
    w3=(322.0 - 13.0*sqrt(70.0))/900.0;
    p1=sqrt( 5.0 - 2.0*sqrt(10.0/7.0) )/3.0;
    p2=sqrt( 5.0 + 2.0*sqrt(10.0/7.0) )/3.0;
    w(1)=w1;   x(1)=0;     y(1)=0.0; z(1)=0.0;
    w(2)=w2;   x(2)=p1;    y(2)=0.0; z(2)=0.0;
    w(3)=w2;   x(3)=-p1;   y(3)=0.0; z(3)=0.0;
    w(4)=w3;   x(4)=p2;    y(4)=0.0; z(4)=0.0;
    w(5)=w3;   x(5)=-p2;   y(5)=0.0; z(5)=0.0;
end

%2 DIMENSIONAL INTEGRATION RULES
%Order 0 polynomial
function [w,x,y,z] = dim_2_order_0
    %Integration scheme integrates linear exactly on unit triangle
    w(1)=1.0/2.0; x(1)=1.0/3.0; y(1)=1.0/3.0; z(1)=0.0;
end

%Order 2 polynomial
function [w,x,y,z] = dim_2_order_2
    %Integration scheme integrate quadratic exactly on unit triangle
    w(1)=1.0/6.0; x(1)=2.0/3.0; y(1)=1.0/6.0; z(1)=0.0;
    w(2)=1.0/6.0; x(2)=1.0/6.0; y(2)=1.0/6.0; z(2)=0.0;
    w(3)=1.0/6.0; x(3)=1.0/6.0; y(3)=2.0/3.0; z(3)=0.0;   
end

%Order 4 polynomial
function [w,x,y,z] = dim_2_order_4
%Integration scheme integrates quartic exactly on unit triangle
    w(1)=0.5*0.2233815896780; x(1)=0.1081030181681; y(1)=0.4459484909160; z(1)=0.0;
    w(2)=0.5*0.2233815896780; x(2)=0.4459484909160; y(2)=0.1081030181681; z(2)=0.0;
    w(3)=0.5*0.2233815896780; x(3)=0.4459484909160; y(3)=0.4459484909160; z(3)=0.0;
    w(4)=0.5*0.1099517436553; x(4)=0.8168475729805; y(4)=0.0915762135098; z(4)=0.0;
    w(5)=0.5*0.1099517436553; x(5)=0.0915762135098; y(5)=0.0915762135098; z(5)=0.0;
    w(6)=0.5*0.1099517436553; x(6)=0.0915762135098; y(6)=0.8168475729805; z(6)=0.0;
end

%Order 5 polynomial
function [w,x,y,z]=dim_2_order_5
    %Integration scheme integrates quintic exactly on unit triangle
    w(1)=0.5*0.1259391805448; x(1)=0.1012865073235; y(1)=0.1012865073235; z(1)=0.0;
    w(2)=0.5*0.1259391805448; x(2)=0.7974269853531; y(2)=0.1012865073235; z(2)=0.0;
    w(3)=0.5*0.1259391805448; x(3)=0.1012865073235; y(3)=0.7974269853531; z(3)=0.0;
    w(4)=0.5*0.1323941527885; x(4)=0.4701420641051; y(4)=0.0597158717898; z(4)=0.0;
    w(5)=0.5*0.1323941527885; x(5)=0.4701420641051; y(5)=0.4701420641051; z(5)=0.0;
    w(6)=0.5*0.1323941527885; x(6)=0.0597158717898; y(6)=0.4701420641051; z(6)=0.0;
    w(7)=0.5*0.225;           x(7)=0.3333333333333; y(7)=0.3333333333333; z(7)=0.0;
end

%Order 7 polynomial
function [w,x,y,z]=dim_2_order_7
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

%3 DIMENSIONAL INTEGRATION RULES
%Order 0 polynomial
function [w,x,y,z] = dim_3_order_0
    %Integration scheme integrates exactly linear on unit tetrahedron
    w(1)=1.0/6.0; x(1)=1.0/4.0; y(1)=1.0/4.0; z(1)=1.0/4.0;
end

%Order 2 polynomial
function [w,x,y,z] = dim_3_order_2
    %Integration scheme integrates exactly quadratic on unit tetrahedron
    a=(5+3*sqrt(5))/20; b=(5-sqrt(5))/20;
    w(1)=1.0/24.0; x(1)=a; y(1)=b; z(1)=b;
    w(2)=1.0/24.0; x(2)=b; y(2)=a; z(2)=b;
    w(3)=1.0/24.0; x(3)=b; y(3)=b; z(3)=a;
    w(4)=1.0/24.0; x(4)=b; y(4)=b; z(4)=b;
end

%Order 4 polynomial
function [w,x,y,z] = dim_3_order_4
    %Integration scheme integrates exactly quartic on unit tetrahedron
    a=0.25;
    w(1)=-0.0131555555555555556;
    x(1)=a; y(1)=a; z(1)=a;
    
    a=0.785714285714285714; b=0.0714285714285714285;
    w(2)=0.00762222222222222222;
    x(2)=a; y(2)=b; z(2)=b;
    w(3)=0.00762222222222222222;
    x(3)=b; y(3)=a; z(3)=b;
    w(4)=0.00762222222222222222;
    x(4)=b; y(4)=b; z(4)=a;
    w(5)=0.00762222222222222222;
    x(5)=b; y(5)=b; z(5)=b;
    
    a=0.399403576166799219; b=0.100596423833200785;
    w(6)=0.0248888888888888889;
    x(6)=a; y(6)=a; z(6)=b;
    w(7)=0.0248888888888888889;
    x(7)=a; y(7)=b; z(7)=b;
    w(8)=0.0248888888888888889;
    x(8)=b; y(8)=b; z(8)=a;
    w(9)=0.0248888888888888889;
    x(9)=b; y(9)=a; z(9)=b;
    w(10)=0.0248888888888888889;
    x(10)=b; y(10)=a; z(10)=a;
    w(11)=0.0248888888888888889;
    x(11)=a; y(11)=b; z(11)=a;
end

end