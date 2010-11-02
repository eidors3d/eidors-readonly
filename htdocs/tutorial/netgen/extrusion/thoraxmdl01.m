subplot(221);

load CT1.mat
img = flipdim(imread('thorax-mdl.jpg'),1); %Keep up direction
colormap(gray(256));
imagesc(img);
 set(gca,'YDir','normal');

thorax = thorax*100;
thorax(:,2) = 512 - thorax(:,2);
rlung = rlung*100;
rlung(:,2) = 512 - rlung(:,2);
llung = llung*100;
llung(:,2) = 512 - llung(:,2);
llung = [llung; 340,375;345,360];

hold on;
plot(thorax(:,1),thorax(:,2),'b','LineWidth',2);
plot(rlung(:,1),rlung(:,2),'b','LineWidth',2);
plot(llung(:,1),llung(:,2),'b','LineWidth',2);
hold off

axis off; axis equal
print_convert thoraxmdl01a.jpg
return

% Shrink the model  for the next step
trunk = trunk*.01;
lung  = lung*.01; lung = flipud(lung(1:3:end,:)); % need counterclockwise shapes
elec_pos = elec_pos*.01;
