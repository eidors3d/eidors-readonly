img_v = img;
img_v.fwd_model.mdl_slice_mapper.npx = 64;
img_v.fwd_model.mdl_slice_mapper.npy = 64;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,0.8];
show_current(img_v, vh.volt(:,1));

axis equal; axis tight; print_convert pig_body04a.jpg

img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];
show_current(img_v, vh.volt(:,1));

axis equal; axis tight; print_convert pig_body04b.jpg
