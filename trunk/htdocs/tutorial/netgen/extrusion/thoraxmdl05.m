img_v.fwd_model.mdl_slice_mapper.npx = 1000;
img_v.fwd_model.mdl_slice_mapper.npy = 1000;
img_v.fwd_model.mdl_slice_mapper.level = PLANE;

% Calculate at high resolution
q = show_current(img_v, vh.volt(:,1));

pic = shape_library('get','adult_male','pic');
imagesc(pic.X, pic.Y, pic.img);
% imgt= flipdim(imread('thorax-mdl.jpg'),1); imagesc(imgt);
colormap(gray(256)); set(gca,'YDir','normal');
hold on

sx = linspace(-.5,.5,15)';
sy = 0.05 + linspace(-.5,.5,15)';
hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',2, 'color','b');
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',2, 'color','b');

axis equal; axis tight; axis off; print_convert thoraxmdl05a.jpg
