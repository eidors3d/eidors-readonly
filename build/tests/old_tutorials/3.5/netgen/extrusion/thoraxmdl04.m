img_v = img;
img_v.fwd_model.mdl_slice_mapper.npx = 64;
img_v.fwd_model.mdl_slice_mapper.npy = 64;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];
show_current(img_v, vh.volt(:,1));

axis tight; axis image; ylim([50,450]); axis off
print_convert thoraxmdl04a.jpg

