img_v.fwd_model.mdl_slice_mapper.npx = 1000;
img_v.fwd_model.mdl_slice_mapper.npy = 1000;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];

% Calculate at high resolution
 q = show_current(img_v, vh.volt(:,1));

% Lower resolution to visualize
img_v.fwd_model.mdl_slice_mapper.npx = 64;
img_v.fwd_model.mdl_slice_mapper.npy = 64;
show_current(img_v, vh.volt(:,1));


sx =-centroid(1) - linspace(-1,1,15)';
sy =-centroid(2) + linspace(-1,1,15)';
hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy);

axis equal; axis tight;  print_convert pig_body05a.jpg
