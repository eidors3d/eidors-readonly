function [ztab,index] = TabulateImageData(z,mask)
% function [ztab,index] = TabulateImageData(z,mask)
% returns a table of time series per row
% index is the mapping from table to (n1,n2) coords and masked state (0
% masked, 1 included)
% mask is optional
% Copyright C. Gomez-Laberge, August 2009
% $Id: $

if nargin == 1
    mask = ones(size(z));
end
[n1,n2,n3] = size(z);
index = [kron(ones(n1,1),(1:n1)'), kron((1:n2)',ones(n2,1)),mask(:)];
% index = [kron((1:n2)',ones(n2,1)),kron(ones(n1,1),(1:n1)'),mask(:)]; % [xy are flipped here]
ztemp = [];
for j = 1:n2
    zj = z(:,j,:);
    ztemp = [ztemp; reshape(zj,[n1,n3])];
end
ztab = ztemp(logical(index(:,3)),:);
end
