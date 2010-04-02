% RPI tank model $Id$

% In 3D, it's important to get the model diameter right, 2D is
imdl.fwd_model.nodes= imdl.fwd_model.nodes*15; % 30 cm diameter

% Estimate the background conductivity
imgs = calc_jacobian_bkgnd(imdl);
vs = fwd_solve(imgs); vs = vs.meas;

pf = polyfit(vh,vs,1);

imdl.jacobian_bkgnd.value = pf(1)*imdl.jacobian_bkgnd.value;
imdl.hyperparameter.value = imdl.hyperparameter.value/pf(1)^2;

subplot(221);
img = inv_solve(imdl, vh, vi);
img.calc_colours.cb_shrink_move = [0.5,0.8,0];
show_fem(img,[1,1]); axis off; axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2  jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps rpi_data04a.png
