function [A,b]=tet_to_inequal(v,e)
% [A,b]=tet_to_inequal(v)
% Given the vertices of a simplex v return a system
% of linear inequalities so that a point x in in 
% the simplex iff Ax <= b
%
% [A,b]=tet_to_inequal(v)
% Given the Nx3 list of vertices v and the Nx4 list of elements e indexing
% v, return a system of linear inequalities so that a point x in in 
% the simplex iff Ax <= b

% (C) 2012 Bill Lionheart. License GPL v2 or v3
% $Id$

if isstr(v) && strcmp(v,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   e = [1 2 3 4];
   ne = 1;
end
ne = size(e,1);
d1 = size(e,2);
d  = size(v,2);
% [d1,d]=size(v(e,:));
if d~=3
    error('Only works for dimension 3 at the moment')
end
if d1 ~= d+1
    error('Simplex needs to have one more vertex than dimension')
end
% the d edges at vertex 1
edges1= v(e(:,2:end)',:) - kron(v(e(:,1),:), ones(3,1));%ones(3,1)*v(e(:,1),:);

dt = my_det(edges1);
idx = dt<0;
if any(idx)
   e(idx, :) = [e(idx,2) e(idx,1) e(idx,3:end)];
   idx3 = reshape(repmat(idx,1,3)',1,[])';
   edges1(idx3,:)= v(e(idx,2:end)',:)-reshape(repmat(v(e(idx,1),:),1,3)',3,[])';
end
idx1 = 1:3*ne; % 1 2 3 4 5 6 ...
idx3 = reshape(circshift(reshape(idx1,3,[]),-1),1,[]); % 2 3 1 5 6 4 ...
idx2 = 1:4*ne; idx2(4:4:end) = [];
A(idx2,:) = -reshape(my_cross(edges1(idx1,:),edges1(idx3,:)),[],3);
b(idx2) = sum(A(idx2,:).*v(e(:,2:end)',:),2);
% One more face not containing vertex 1
edges3= [v(e(:,3)',:)-v(e(:,2)',:);v(e(:,4)',:)-v(e(:,2)',:)];
idx4 = 4:4:4*ne;
A(idx4,:) = reshape(my_cross(edges3(1:end/2,:),edges3(end/2+1:end,:)),[],3);
b(idx4) = sum(A(idx4,:).*v(e(:,2)',:),2);
b=b';
end

function c = my_cross(a,b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) 
     a(:,3).*b(:,1)-a(:,1).*b(:,3)
     a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end  

function d = my_det(a)
ln = size(a,1);
c1 = 1:3:ln;
c2 = 2:3:ln;
c3 = 3:3:ln;
d = a(c1,1).*a(c2,2).*a(c3,3) + ...
    a(c2,1).*a(c3,2).*a(c1,3) + ...
    a(c3,1).*a(c1,2).*a(c2,3) - ...
    a(c3,1).*a(c2,2).*a(c1,3) - ...
    a(c2,1).*a(c1,2).*a(c3,3) - ...
    a(c1,1).*a(c3,2).*a(c2,3);
end

function do_unit_test
simple_test
test_3d
end

function test_3d
imdl = mk_common_model('a3cr',16);
fmdl = imdl.fwd_model;
fmdl = linear_reorder(fmdl,1);
[A1,b1] = tet_to_inequal(fmdl.nodes,fmdl.elems);
v = fmdl.nodes(fmdl.elems(1,:),:);
[A2,b2] = tet_to_inequal(v);
unit_test_cmp('Parallel vs serial', A1((1:4),:),A2)
unit_test_cmp('Parallel vs serial', b1(1:4),b2)
v = fmdl.nodes(fmdl.elems(2,:),:);
[A2,b2] = tet_to_inequal(v);
unit_test_cmp('Parallel vs serial', A1(4*(2-1)+(1:4),:),A2)
unit_test_cmp('Parallel vs serial', b1(4*(2-1)+(1:4)),b2)
v = fmdl.nodes(fmdl.elems(57,:),:);
[A2,b2] = tet_to_inequal(v);
unit_test_cmp('Parallel vs serial', A1(4*(57-1)+(1:4),:),A2)
unit_test_cmp('Parallel vs serial', b1(4*(57-1)+(1:4)),b2)
end

function simple_test
v1=[0,0,0;eye(3)];
[A1,b1] =  tet_to_inequal(v1);
out = all( A1*[0.1;0.1;0.1]-b1<0);
correct = 1;
unit_test_cmp('Point in unit tetrahedron', out, correct)
end