% Show moving ball 
% $Id: moving_ball_16x1_03.m,v 1.2 2007-09-25 11:21:43 aadler Exp $
show_fem(ng_mdl_16x1_coarse);
crop_model(gca, inline('z>9','x','y','z'))
view(0,70)
[xs,ys,zs] = sphere(n);

surf(xs,ys,zs);
print -r100 -dpng build_single_plane03a.png
