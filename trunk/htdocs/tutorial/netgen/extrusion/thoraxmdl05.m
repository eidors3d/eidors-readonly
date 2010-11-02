img_v.fwd_model.mdl_slice_mapper.npx = 200;
img_v.fwd_model.mdl_slice_mapper.npy = 200;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];

% Calculate at high resolution
 q = show_current(img_v, vh.volt(:,1));

img = flipdim(imread('thorax-mdl.jpg'),1); imagesc(img);
colormap(gray(256)); set(gca,'YDir','normal');


hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',2);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',2);

axis equal; axis tight; axis off; print_convert pig_body06a.jpg



