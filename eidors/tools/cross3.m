function c = cross3(a,b)
%CROSS3  3D cross parallel cross product
%  C = CROSS3(A,B) calculates the cross product between rows of matrices A
%  and B, both of which must be Nx3. C is also Nx3 such that 
%    C(i,:) = cross(A(i,:),B(i,:))
% but it is calculated faster.

% (C) 2013 Bartlomiej Grychtol. License: GPL version 2 or 3.
% $Id$

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
     a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
     a(:,1).*b(:,2)-a(:,2).*b(:,1)];