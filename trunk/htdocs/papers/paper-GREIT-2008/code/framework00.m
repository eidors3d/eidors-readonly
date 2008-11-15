% fwd_model $Id$

load ng_mdl_16x1_fine

subplot(121)
show_fem(ng_mdl_16x1_fine);
view([-5 28]);
crop_model(gca, inline('x+2*z>20','x','y','z'))

subplot(122)
show_fem(ng_mdl_16x1_fine);
view(0,0);
crop_model(gca, inline('y>-10','x','y','z'))
set(gca,'Xlim',[-4,4],'Zlim',[-2,2]+5);

print -dpng -r150 framework00a.png
