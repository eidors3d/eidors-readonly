% Show moving ball 
% $Id: moving_ball_16x1_03.m,v 1.1 2007-09-24 20:53:24 aadler Exp $
show_fem(ng_mdl_16x1_coarse);
crop_model(gca, inline('z>9','x','y','z'))
view(0,70)
print -r100 -dpng build_single_plane02a.png
