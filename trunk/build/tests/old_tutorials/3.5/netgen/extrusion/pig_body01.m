subplot(221);

load CT2.mat
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

% Shrink the model  for the next step
trunk = trunk*.01;
lung  = lung*.01; lung = flipud(lung(1:3:end,:)); % need counterclockwise shapes
elec_pos = elec_pos*.01;
