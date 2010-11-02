load CT1.mat
img = flipdim(imread('thorax-mdl.jpg'),1); %Keep up direction
imagesc(img);
colormap(gray(256)); set(gca,'YDir','normal');

hold on;
plot(thorax(:,1),thorax(:,2),'b','LineWidth',2);
plot(rlung(:,1),rlung(:,2),'b','LineWidth',2);
plot(llung(:,1),llung(:,2),'b','LineWidth',2);
hold off

subplot(221); axis off; axis equal
print_convert thoraxmdl01a.jpg
