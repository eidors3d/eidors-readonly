img_v.fwd_model.mdl_slice_mapper.npx = 128;
img_v.fwd_model.mdl_slice_mapper.npy = 128;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];

 q = show_current(img_v, vh.volt(:,1));

sx =-centroid(1) - linspace(-1,1,15)';
sy =-centroid(2) + linspace(-1,1,15)';
hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy);

axis tight; print_convert pig_body05a.jpg
