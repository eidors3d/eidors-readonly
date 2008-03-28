% Simulate Moving Ball $Id: square_mesh04.m,v 1.1 2008-03-28 16:18:50 aadler Exp $

% get ng_mdl_16x1_vfine from data_contrib section of web page
n_sims= 20;
 load ng_mdl_16x1_vfine.mat; fmdl= ng_mdl_16x1_vfine;
%load ng_mdl_16x1_coarse.mat; fmdl= ng_mdl_16x1_coarse;
[vh,vi,xyzr_pt]= simulate_3d_movement( n_sims, fmdl);

clf;
show_fem(fmdl)
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
print -r100 -dpng square_mesh04a.png;
view(-12,4)
print -r100 -dpng square_mesh04b.png;

