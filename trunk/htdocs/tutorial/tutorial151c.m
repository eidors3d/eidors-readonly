% Difference imaging result
% $Id: tutorial151c.m,v 1.3 2006-11-17 03:57:30 aadler Exp $

% imdl is loaded from file tutorial151_model.mat
imdl.reconst_type= 'difference';
imdl.solve=        @np_inv_solve;
imdl.jacobian_bkgnd.value= backgnd;
imdl.hyperparameter.value= 1e-2;
img_diff= inv_solve(imdl, v_homg, v_targ);

clf;
axes('position',[.1,.1,.65,.8]);
   show_fem(img_diff); view(-33,20);
axes('position',[.8,.1,.15,.8]);
   show_slices(img_diff,5);

print -r75 -dpng tutorial151c.png;
