function show_polygonsValues(imdl,data,colors,edgeColor)

if nargin<3
    colors= 'jet';
    edgeColor= 'none';
elseif nargin<4
    edgeColor= 'none';
end

polygonsX= imdl.inv_solve.polygonsX;
polygonsY= imdl.inv_solve.polygonsY;

if size(data,1)>1 % data must be a row vector
    data= data';
end
if size(polygonsX,2)~= size(data,2)
    polygonsX= polygonsX';
end
if size(polygonsY,2)~= size(data,2)
    polygonsY= polygonsY';
end

fill(polygonsX,polygonsY,data,'EdgeColor',edgeColor);
hold on;
% colorbar
if ischar(colors)
    colormap(colors)
else
    set(gcf,'Colormap',colors)
end


end