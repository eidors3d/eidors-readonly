% Simulate Moving Ball $Id$

n_sims= 20;
f_mdl = mk_library_model('cylinder_16x1el_vfine');
f_mdl.stimulation = mk_stim_patterns(16,1,'{ad}','{ad}',{},1);
[vh,vi,xyzr_pt]= simulate_3d_movement( n_sims, f_mdl);

clf;
show_fem(f_mdl)
crop_model(gca, inline('x-z<-8','x','y','z'))

hold on
[xs,ys,zs]=sphere(10);
for i=1:n_sims
   xp=xyzr_pt(1,i); yp=xyzr_pt(2,i);
   zp=xyzr_pt(3,i); rp=xyzr_pt(4,i);
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',[.4,0,.4],'FaceColor',[.2,0,.2]);
end
zofs=1/3;
plot3(xpts*15,ypts*15,(zpts+zofs)*15,'b');
hold off

axis equal
view(-23,44)
print_convert square_mesh04a.png '-density 75'
view(-12,4)
print_convert square_mesh04b.png '-density 75'

