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

% plot(imdl.fwd_model.misc.elec_posn(:,2),imdl.fwd_model.misc.elec_posn(:,3),'*k','linewidth',2,'MarkerSize',10)
% xlabel('X (m)','fontsize',20,'fontname','Times');
% ylabel('Y (m)','fontsize',20,'fontname','Times');
% axis equal; axis tight;
% set(gca,'fontsize',20,'fontname','Times')


end