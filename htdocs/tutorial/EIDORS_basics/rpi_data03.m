% RPI tank model $Id$

% simple inverse model -> replace fields to match this model
imdl.fwd_model = fmdl; 

subplot(221);
img = inv_solve(imdl , vh, vi);
img.calc_colours.cb_shrink_move = [0.5,0.8,0];
show_fem(img,[1,1]); axis off; axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2  jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps rpi_data03a.png
