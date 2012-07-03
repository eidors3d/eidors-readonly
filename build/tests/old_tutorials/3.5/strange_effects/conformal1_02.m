% Voltage distribution $Id$

img.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img);

imgv = rmfield(img,'elem_data');
imgv.node_data = vv.volt;

show_fem( imgv ); axis image

%print -dpng -r125 rpi_data01a.png
print -depsc2 jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps conformal1_02a.png

img2.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img2);

imgv = rmfield(img2,'elem_data');
imgv.node_data = vv.volt;

show_fem( imgv ); axis image
%print -dpng -r125 rpi_data01a.png
print -depsc2 jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps conformal1_02b.png
