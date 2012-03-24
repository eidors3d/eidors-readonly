imgc= img;
imgc.fwd_model.mdl_slice_mapper.npx = 128;
imgc.fwd_model.mdl_slice_mapper.npy = 200;
imgc.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
hh=show_fem(imgc);
set(hh,'EdgeColor',[1,1,1]*.75);

 hold on;
q = show_current(imgc,vh.volt);
quiver(q.xp,q.yp, q.xc,q.yc,15,'b','LineWidth',1);
axis([-.2,.2,0.8,1.05]);
hold off;

title(sprintf('current near electrode:  zc = %5.3f',fmdl.electrode(1).z_contact));
print_convert contact_impedance03a.png
