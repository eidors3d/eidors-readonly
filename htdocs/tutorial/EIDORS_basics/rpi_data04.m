% RPI tank model $Id$

% In 3D, it's important to get the model diameter right, 2D is
imdl.fwd_model.nodes= imdl.fwd_model.nodes*15; % 30 cm diameter

% Estimate the background conductivity
imgs = mk_image(imdl);
vs = fwd_solve(imgs); vs = vs.meas;

pf = polyfit(vh,vs,1);

imdl.jacobian_bkgnd.value = pf(1)*imdl.jacobian_bkgnd.value;
imdl.hyperparameter.value = imdl.hyperparameter.value/pf(1)^2;

img = inv_solve(imdl, vh, vi);
img.calc_colours.cb_shrink_move = [0.5,0.8,0];
show_fem(img,[1,1]); axis off; axis image

print_convert('rpi_data04a.png','-density 60');

