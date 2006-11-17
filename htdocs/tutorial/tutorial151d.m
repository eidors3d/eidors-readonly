% Difference imaging result
% $Id: tutorial151d.m,v 1.2 2006-11-17 03:57:30 aadler Exp $

% imdl is loaded from file tutorial151_model.mat
imdl.reconst_type= 'static';
imdl.solve=        @tutorial151_nonlinearGN;
imdl.parameters.term_tolerance= 1e-5;
imdl.parameters.max_iterations= 10;
% special parameter for this model
imdl.tutorial151_nonlinearGN.init_backgnd= backgnd;

imdl.hyperparameter.value= 3e-1;
img_diff= inv_solve(imdl, v_targ);

clf;
axes('position',[.1,.1,.65,.8]);
   show_fem(img_diff); view(-33,20);
axes('position',[.8,.1,.15,.8]);
   show_slices(img_diff,5);

print -r75 -dpng tutorial151d.png;
