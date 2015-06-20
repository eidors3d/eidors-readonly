% positions of where to plot
yxposns = [30 30 15;
           22 50 40]*2;

% Show image
clf; axes('position',[0.05,0.5,0.25,0.45]);
img1= imgr; img1.elem_data = imgr.elem_data(:,1);
img1.show_slices.img_cols= 1;
show_slices(img1);
hold all;
for i = 1:size(yxposns,2);
    plot(yxposns(2,i),yxposns(1,i),'s','LineWidth',2);
end
hh=text(62,14,'V','FontSize',10,'FontWeight','bold');
hh=text(62,115,'D','FontSize',10,'FontWeight','bold');
hh=text( 2,64,'R','FontSize',10,'FontWeight','bold');
hh=text(120,64,'L','FontSize',10,'FontWeight','bold');
hold off;

% Show plots
imgs =-calc_slices(imgall);
axes('position',[0.37,0.62,0.33,0.28]);

imgs = permute(imgs,[3,1,2]);
taxis =  (0:size(imgs,1)-1)/13; % frame rate = 13
hold all
for i = 1:size(yxposns,2);
    plot(taxis,imgs(:,yxposns(1,i),yxposns(2,i)),'LineWidth',2);
end
hold off
set(gca,'ytick',[]);
axis tight
xlim([0 7]);
ylabel '\Delta Resistivity ->'
xlabel 'Time (s)'

 print_convert neonate_ex04a.png
