function index = Coord2Index(coord,dims)
% NOTE COORDS SHOWN AS [X,Y] FROM MATLAB IMAGE ARE ACTUALLY [Y,X] AND SO
% ARE TRANSPOSED HERE.
% coord = [x1,y1;x2,y2;...;xm,ym];
% dims = [number of x pixels,number of y pixels];
% Copyright C. Gomez-Laberge, August 2009
% $Id$
m = size(coord,1);
index = zeros(m,1);
xdim = dims(1);
for i = 1:m
    index(i) = coord(i,2)+xdim*(coord(i,1)-1);
end
