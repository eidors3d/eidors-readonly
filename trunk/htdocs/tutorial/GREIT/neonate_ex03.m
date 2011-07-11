vh = mean(vv,2);        % reference is average
img = inv_solve(imdl,vh,vv);
img.calc_colours.ref_level=0;
img.calc_colours.npoints  =32;

% Yposns of where to plot
yposns = [10:5:25];

% Show image
clf; axes('position',[0.05,0.5,0.25,0.45]);
img1= img; img1.elem_data = img.elem_data(:,45);
show_slices(img1);
hold on;
plot(10,yposns,'s','LineWidth',5);
hold off;

% Show plots
imgs = calc_slices(img);
axes('position',[0.30,0.6,0.65,0.25]);

imgs = permute(imgs,[3,1,2]);
plot(imgs(:,yposns,10),'LineWidth',2);
axis tight

print_convert neonate_ex03a.png
