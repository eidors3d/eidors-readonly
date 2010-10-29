subplot(221);

load ~/docs/eidors/htdocs/tutorial/netgen/extrusion/CT2.mat
img = imread('pig-thorax.jpg');
colormap(gray(256));
imagesc(-67+[1:371],25+[1:371],img);

hold on;
plot(372-trunk(:,1),trunk(:,2),'m*-');
plot(372-lung(:,1),lung(:,2),'r*-');
hh=plot(372-elec_pos(:,1),elec_pos(:,2), 'b.'); set(hh,'MarkerSize',20);
hold off

axis off; axis equal
print_convert pig_body01.jpg
