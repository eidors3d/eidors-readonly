% Difference imaging result
% $Id$

% imdl is loaded from file tutorial151_model.mat
imdl.reconst_type= 'difference';
imdl.solve=        @inv_solve_diff_GN_one_step;
imdl.jacobian_bkgnd.value= backgnd;
imdl.hyperparameter.value= 1e-2;
img_diff= inv_solve(imdl, v_homg, v_targ);

clf;
axes('position',[.1,.1,.65,.8]);
   show_fem(img_diff); view(-33,20);
axes('position',[.8,.1,.15,.8]);
   show_slices(img_diff,5);

print_convert tutorial151c.png '-density 75'
