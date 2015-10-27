%%BUG. NEED TO FIND IMGALL.

% positions of where to plot
yposns = [20 20 45 45];
xposns = [20 45 20 45];

% Show image
clf; axes('position',[0.05,0.5,0.25,0.45]);
img1= imgall; img1.elem_data = imgall.elem_data(:,45);
show_slices(img1);
hold all;
for i = 1:4
    plot(xposns(i),yposns(i),'s','LineWidth',5);
end
hold off;

% Show plots
imgs = calc_slices(imgall);
axes('position',[0.32,0.6,0.63,0.28]);

imgs = permute(imgs,[3,1,2]);
taxis =  (0:size(imgs,1)-1)/13; % frame rate = 13
hold all
for i = 1:4
    plot(taxis,imgs(:,yposns(i),xposns(i)),'LineWidth',2);
end
hold off
set(gca,'ytick',[]);
xlim([0 16]);

print_convert neonate_ex03a.png
