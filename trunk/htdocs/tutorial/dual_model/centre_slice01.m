% Simulate Moving Ball - Helix $Id: centre_slice01.m,v 1.5 2008-03-25 19:21:49 aadler Exp $

% get ng_mdl_16x2_vfine from data_contrib section of web page
n_sims= 20;
load ng_mdl_16x2_vfine.mat; fmdl= ng_mdl_16x2_vfine;
%imdl= mk_common_model('b3cr',[16,2]); fmdl= imdl.fwd_model;
[vh,vi,xyzr_pt]= simulate_3d_movement( n_sims, fmdl);

show_fem(fmdl)
crop_model(gca, inline('x-z<-15','x','y','z'))
%crop_model(gca, inline('x-z<-0.5','x','y','z'))
view(-23,14)

hold on
[xs,ys,zs]=sphere(10);
for i=1:n_sims
   xp=xyzr_pt(1,i); yp=xyzr_pt(2,i);
   zp=xyzr_pt(3,i); rp=xyzr_pt(4,i);
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',[1,0,1],'FaceColor',[1,0,1]);
end
hold off

print -r100 -dpng centre_slice01a.png;
