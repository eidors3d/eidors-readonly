% Cheap static solver $Id$

% Do a diff solve with respect to simulated data
imgs = mk_image(imdl);
vs = fwd_solve(imgs); vs = vs.meas;


imdl.hyperparameter.value = 1.00;
img = inv_solve(imdl, vs, vi);
img.elem_data = pf(1) + img.elem_data*0.5;
img.calc_colours.cb_shrink_move = [0.5,0.8,0.05];
show_fem(img,[1,1]); axis off; axis image

print_convert('rpi_data05a.png','-density 60');
