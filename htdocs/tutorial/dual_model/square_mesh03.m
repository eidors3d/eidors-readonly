% Dual models $Id$

load ng_mdl_16x1_coarse; f_mdl= ng_mdl_16x1_coarse;

subplot(122)
show_fem(f_mdl);  % fine model
crop_model(gca, inline('x-z<-8','x','y','z'))

% Map coarse model geometry
zofs=1/3;
hold on
plot3(xpts*15,ypts*15,(zpts+zofs)*15,'b');
hold off

axis(15*[-1.1,+1.1,-1.1,+1.1,zofs-0.4,zofs+0.4]);
view(-47,28); axis square

print -r150 -dpng square_mesh03a.png
