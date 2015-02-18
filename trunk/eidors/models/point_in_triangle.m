function out = point_in_triangle(P,E,V,epsilon, str)
%POINT_IN_TRIANGLE tests points for membership in triangles
% POINT_IN_TRIANGLE(P, E, V)
% POINT_IN_TRIANGLE(P, E, V, epsilon)
% Inputs:
%   P       - [p x 2] or [p x 3] list of point coordinates
%   V       - [v x 2] or [v x 3] list of triangle vertices
%   E       - [e x 3] matrix specifying the vertices V consituting each of the
%               e triangles
%   epsilon - threshold to counteract numerical instability:
%               epsilon > 0 makes the triangle bigger so points on edges
%                   and vertices are correctly classified as inside. May
%                   contain false positives.
%               epsilon = 0 may miss some points on edges and vertices due 
%                   to numerical issues.
%               epsilon < 0 makes the trianle smaller so points on edges
%                   and vertices are classified as outside. May miss more
%                   points than intended.
%               Default: epsilon = eps
%
% POINT_IN_TRIANGLE(..., 'match') specifies that each point is to be tested
%   only against its corresponding triangle (number of points must equal
%   number of triangles). Option has no effect in 2D.


% (C) 2012-2015 Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

switch nargin
    case {1 2 3}
        epsilon = eps;
        str = [];
    case 4
        if ischar(epsilon)
            str = epsilon;
            epsilon = eps;
        else
            str = [];
        end
end


match = strcmp(str,'match');

switch size(P,2)
    case 2
        [u, v] = point_in_triangle_2d(P,E,V);
    case 3
        [u, v] = point_in_triangle_3d_wrapper(P,E,V,match);
        
    otherwise
        error('EIDORS:WrongInput','Points must be 2D or 3D');
end

out = u >= -epsilon & v >= -epsilon & (u+v-epsilon) <= 1; 

function [u, v] = point_in_triangle_3d_wrapper(P,E,V, match)
    nPts = size(P,1);
    nTri = size(E,1);
    
    if nPts == 1 || match
        [u, v] = point_in_triangle_3d(P,E,V);
    else
        u = zeros(nPts, nTri);
        v = zeros(nPts, nTri);
        for i = 1:nPts
            [u(i,:), v(i,:)] = point_in_triangle_3d(P(i,:),E,V);
        end

    end

%-------------------------------------------------------------------------%
% Decide if point P is in triangles E indexing nodes V in 3D
function [u, v] = point_in_triangle_3d(p, E, V)
%http://www.blackpawn.com/texts/pointinpoly/default.html
% vectors
v0 = V(E(:,3),:) - V(E(:,1),:);
v1 = V(E(:,2),:) - V(E(:,1),:);
v2 = bsxfun(@minus,p, V(E(:,1),:));

% dot products
dot00 = dot(v0, v0, 2);
dot01 = dot(v0, v1, 2);
dot02 = dot(v0, v2, 2);
dot11 = dot(v1, v1, 2);
dot12 = dot(v1, v2, 2);

% barycentric coordinates
invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;


%-------------------------------------------------------------------------%
% Decide if points P are in triangles E indexing nodes V in 2D
function [u, v] = point_in_triangle_2d(P,E,V)
    X = reshape(V(E,1),size(E));
    Y = reshape(V(E,2),size(E));
    T = [ bsxfun(@minus, X(:,1:2), X(:,3)) bsxfun(@minus, Y(:,1:2), Y(:,3))];
    invdetT = 1./(T(:,1).*T(:,4) - T(:,2).*T(:,3));
    nPts = size(P,1);
    u = ( repmat(Y(:,2)-Y(:,3),1,nPts).*bsxfun(@minus,P(:,1)',X(:,3)) ...
        + repmat(X(:,3)-X(:,2),1,nPts).*bsxfun(@minus,P(:,2)',Y(:,3)) )...
        .* repmat(invdetT,1,nPts);
    v = ( repmat(Y(:,3)-Y(:,1),1,nPts).*bsxfun(@minus,P(:,1)',X(:,3)) ...
        + repmat(X(:,1)-X(:,3),1,nPts).*bsxfun(@minus,P(:,2)',Y(:,3)) )...
        .* repmat(invdetT,1,nPts);
    u = u';
    v = v';
    

