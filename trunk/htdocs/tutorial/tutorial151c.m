% Difference imaging result
% $Id: tutorial151c.m,v 1.1 2006-11-15 17:16:27 aadler Exp $

% imdl is loaded from file tutorial151_model.mat
imdl.reconst_type= 'difference';
imdl.solve=        @np_inv_solve;
imdl.jacobian_bkgnd.value= backgnd;
img_diff= inv_solve(imdl, v_homg, v_targ);

clf;
axes('position',[.1,.1,.65,.8]);
   show_fem(sim_img); view(-33,20);
axes('position',[.8,.1,.15,.8]);
   show_slices(sim_img,5);

print -r75 -dpng tutorial151a.png;
