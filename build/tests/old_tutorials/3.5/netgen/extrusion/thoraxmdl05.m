img_v.fwd_model.mdl_slice_mapper.npx = 1000;
img_v.fwd_model.mdl_slice_mapper.npy = 1000;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];

% Calculate at high resolution
 q = show_current(img_v, vh.volt(:,1));

imgt= flipdim(imread('thorax-mdl.jpg'),1); imagesc(imgt);
colormap(gray(256)); set(gca,'YDir','normal');

sx = 250 - linspace(-130,130,15)';
sy = 250 + linspace(-130,130,15)';

hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',2);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',2);

axis equal; axis tight; axis off; print_convert thoraxmdl05a.jpg
