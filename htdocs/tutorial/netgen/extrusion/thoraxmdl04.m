img_v = img;
img_v.fwd_model.mdl_slice_mapper.npx = 64;
img_v.fwd_model.mdl_slice_mapper.npy = 64;
img_v.fwd_model.mdl_slice_mapper.level = PLANE;
q = show_current(img_v, vh.volt(:,1));
quiver(q.xp,q.yp, q.xc,q.yc,10,'b');
axis tight; axis image; ylim([-1 1]);axis off
print_convert thoraxmdl04a.jpg

