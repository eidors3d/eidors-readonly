function mdlSurf= show_femTopography(fmdl)
% load('matFiles/fmdlTournemireTopo','fmdl');

En = fmdl.nodes(:,1) == max(fmdl.nodes(:,1));
E  = all(En(fmdl.boundary),2);

Wn = fmdl.nodes(:,1) == min(fmdl.nodes(:,1));
W  = all(Wn(fmdl.boundary),2);


Sn = fmdl.nodes(:,2) == min(fmdl.nodes(:,2));
S  = all(Sn(fmdl.boundary),2);

Nn = fmdl.nodes(:,2) == max(fmdl.nodes(:,2));
N  = all(Nn(fmdl.boundary),2);


Fn = fmdl.nodes(:,3) == min(fmdl.nodes(:,3));
F  = all(Fn(fmdl.boundary),2);

Tn = (fmdl.nodes(:,3) > min(fmdl.nodes(:,3))) & (fmdl.nodes(:,3) < 536);
T  = all(Tn(fmdl.boundary),2);


walls = F | N | E | W | S | T;

tmp.nodes = fmdl.nodes;
tmp.type  = 'fwd_model';
tmp.elems = fmdl.boundary(~walls ,:);
tmp.boundary = tmp.elems;
tmp.electrode = fmdl.electrode;
h1 = show_fem(tmp);
ch = get(gca,'Children');
for i = 1:length(ch)-1
      set(ch(i),'MarkerFaceColor',get(ch(i),'MarkerEdgeColor'),...
         'MarkerSize',2,'LineStyle','none');
end


mdlSurf= tmp; 
mdlSurf = cleanup_nodes(mdlSurf);
 
set(h1,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.8) %,'FaceAlpha',0.7
set(h1,'EdgeLighting','none');
set(h1,'LineWidth',0.1);

hold on
% walls = W | S;
tmp.elems = fmdl.boundary(walls ,:);
tmp.boundary = tmp.elems;
% tmp = rmfield(tmp,'electrode');
h2 = show_fem(tmp);
set(h2,'FaceLighting','flat','FaceColor','flat','AmbientStrength',0.8) % ,'FaceAlpha',0.7
set(h2,'EdgeLighting','none');
set(h2,'LineWidth',0.1);
set(h2,'FaceColor',[.5 .5 .5]);
hold off

light('Position',[-1 0 0],'Style','infinite');

% caxis([650 850])
xlabel('X (m)','fontsize',16,'fontname','Times')
ylabel('Y (m)','fontsize',16,'fontname','Times')
zlabel('Z (m)','fontsize',16,'fontname','Times')
set(gca,'fontsize',16,'fontname','Times')


return
% print an oversize png
print_convert('tmp.png','-density 1200')
print_convert('tmp.png','-density 1000')

% downsample (introduces anti-aliasing)
I = imread('tmp.png');
I = imresize(I, 1/4, 'bilinear');
imwrite(I,'FigureFEMTournemireTopo.png','png')
delete('tmp.png');

end


function mdl = cleanup_nodes(mdl)
used_nodes = unique(mdl.elems);
un = zeros(1,length(mdl.nodes));
un(used_nodes) = 1:length(used_nodes);
mdl.elems = un(mdl.elems);
mdl.nodes = mdl.nodes(used_nodes,:);
end
