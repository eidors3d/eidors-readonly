img_v = img;
% Stimulate between elecs 16 and 5 to get more interesting pattern
img_v.fwd_model.stimulation(1).stim_pattern = sparse([16;5],1,[1,-1],16,1);
img_v.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img_v);

img_v = rmfield(img, 'elem_data');
img_v.node_data = vh.volt(:,1);
img_v.calc_colours.npoints = 256;
imgs = calc_slices(img_v,PLANE);


clf
imagesc(pic.X, pic.Y, pic.img); colormap(gray(256)); set(gca,'YDir','normal');
hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',2);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',2);

[x y] = meshgrid( linspace(pic.X(1), pic.X(2),size(imgs,1)), ...
                  linspace(pic.Y(2), pic.Y(1),size(imgs,2)));
hold on;
contour(x,y,imgs,31);
hh= findobj('Type','patch'); set(hh,'LineWidth',2)

hold off; axis off; axis equal; %ylim([50,450]);
print_convert thoraxmdl06a.jpg
