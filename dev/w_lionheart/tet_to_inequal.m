function [A,b]=tet_to_inequal(v)
% [A,b]=tet_to_inequal(v)
% Given the vertices of a simplex v return a system
% of linear inequalities so that a point x in in 
% the simplex iff Ax <= b

% (C) 2012 Bill Lionheart. License GPL v2 or v3
% $Id$

if isstr(v1) && strcmp(v1,'UNIT_TEST'); do_unit_test; return; end

[d1,d]=size(v);
if d~=3
    error('Only works for dimension 3 at the moment')
end
if d1 ~= d+1
    error('Simplex needs to have one more vertex than dimension')
end
% the d edges at vertex 1
edges1= v(2:end,:)-ones(3,1)*v(1,:);
if det(edges1)<0
    %swap two vertices if neg oriented
    v = [v(2,:);v(1,:);v(3:end,:)];
    edges1= v(2:end,:)-ones(3,1)*v(1,:);
end
for i=1:d
    A(i,:) = -cross(edges1(i,:),edges1( mod(i,d)+1,:));
    b(i) = A(i,:)*v(i+1,:)';
end
% One more face not containing vertex 1
edges2= [v(3,:)-v(2,:);v(4,:)-v(2,:)];
A(4,:) = cross(edges2(1,:),edges2(2,:));
b(4) = A(4,:)*v(2,:)';
b=b';

function do_unit_test 
v1=[0,0,0;eye(3)];
[A1,b1] =  tet_to_inequal(v1);
out = all( A1*[0.1;0.1;0.1]-b1<0)
correct = 1;
unit_test_cmp('test1', out, correct)
end