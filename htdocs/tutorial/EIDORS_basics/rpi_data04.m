% Absolute reconstructions $Id$

imdl.solve = @GN_abs_solver;

subplot(221);
imdl.parameters.max_iteration = 1;
img = inv_solve(imdl , vi);
show_fem(img,[0,1]); axis off; axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2  jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps rpi_data03a.png
