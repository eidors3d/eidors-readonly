% Solve 2D and 3D model $Id: two_and_half_d03.m,v 1.1 2008-03-19 00:04:28 aadler Exp $

c2f= mk_coarse_fine_mapping( f_mdl, c_mdl );

imdl.fwd_model.coarse2fine = c2f;
img2= inv_solve(imdl, vh, vi);
subplot(133)
show_fem(img2); view(-62,28)

% Original
subplot(131)
show_fem(demo_img); view(-62,28)

print -r125 -dpng two_and_half_d02a.png
