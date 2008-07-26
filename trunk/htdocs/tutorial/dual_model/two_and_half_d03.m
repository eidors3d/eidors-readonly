% Solve 2D and 3D model $Id$

c2f= mk_coarse_fine_mapping( f_mdl, c_mdl );

imdl.fwd_model.coarse2fine = c2f;
img2= inv_solve(imdl, vh, vi);
img2.elem_data= c2f*img2.elem_data;
subplot(143)
show_fem(img2); view(-62,28)

% 2.5D reconstruct onto coarse model
subplot(144)
img3= inv_solve(imdl, vh, vi);
img3.fwd_model= c_mdl;
show_fem(img3);
set(gca,'ZLim',[0,3]); axis equal; view(-62,28)

print -r150 -dpng two_and_half_d03a.png
