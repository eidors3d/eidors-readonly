% Simulate data $Id: dual_model03.m,v 1.1 2008-03-16 01:33:12 aadler Exp $

img1= inv_solve(imdl(1), vh, vi);
subplot(221)
show_fem(img1);

img2= inv_solve(imdl(2), vh, vi);
c2f= img2.fwd_model.coarse2fine;
img2.elem_data= c2f * img2.elem_data;
subplot(222)
show_fem(img2);

print -r125 -dpng dual_model03a.png;
