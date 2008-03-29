% Reconstruct images $Id: simulate_move2_06.m,v 1.1 2008-03-29 15:48:11 aadler Exp $

imdl= mk_common_model('c2c2',16);
img= inv_solve(imdl, vh, vi);
subplot(121); show_fem(img); axis image
subplot(122); show_slices(img)

print -dpng -r125 simulate_move2_06a.png
