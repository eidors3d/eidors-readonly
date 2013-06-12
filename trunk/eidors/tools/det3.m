function D = det3(a)
%DET3  Parallel determinant for 3x3 matrices
%  D = DET3(A) calculates the determinant for each 3 rows of A. A must be 
%  3N x 3. D is N x 1, such that
%    D(i) = det( A( 3*(i-1)+(1:3) , 1:3 ) )
% but it is calculated faster.

% (C) 2013 Bartlomiej Grychtol. License: GPL version 2 or 3.
% $Id$
ln = size(a,1);
c1 = 1:3:ln;
c2 = 2:3:ln;
c3 = 3:3:ln;
D = a(c1,1).*a(c2,2).*a(c3,3) + ...
    a(c2,1).*a(c3,2).*a(c1,3) + ...
    a(c3,1).*a(c1,2).*a(c2,3) - ...
    a(c3,1).*a(c2,2).*a(c1,3) - ...
    a(c2,1).*a(c1,2).*a(c3,3) - ...
    a(c1,1).*a(c3,2).*a(c2,3);
