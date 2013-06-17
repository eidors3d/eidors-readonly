% Deformation $Id$

xy = interp_mesh(img2.fwd_model);

% Create anisotropic conductivity.
clear conduct
   conduct(:,1,2,2) = abs(xy(:,1));
   conduct(:,1,1,1) = 1;
   conduct(:,1,1,2) = 0;
   conduct(:,1,2,1) = 0;
img2.elem_data = conduct;


img2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img2);

img_v.node_data = vh.volt;
hh=show_fem(img_v);
set(hh,'EdgeColor',[1,1,1]*.75);

q=  show_current(img2,vh.volt);
hold on;
sy =  linspace(1.2,6.8,10); sx =  0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;

axis image
print_convert anisotropy1_03a.png '-density 125'
