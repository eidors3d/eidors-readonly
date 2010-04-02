% Cheap static solver $Id$

% Do a diff solve with respect to simulated data
imgs = calc_jacobian_bkgnd(imdl);
vs = fwd_solve(imgs); vs = vs.meas;


imdl.hyperparameter.value = 18;
img = inv_solve(imdl, vs, vi);
img.elem_data = pf(1) + img.elem_data*0.5;
img.calc_colours.cb_shrink_move = [0.5,0.8,0];
show_fem(img,[1,1]); axis off; axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2  jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps rpi_data05a.png
