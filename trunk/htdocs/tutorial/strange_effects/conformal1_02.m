% Voltage distribution $Id$

img.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img);

imgv = rmfield(img,'elem_data');
imgv.node_data = vv.volt;

show_fem( imgv ); axis image

print_convert conformal1_02a.png '-density 125'

img2.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img2);

imgv = rmfield(img2,'elem_data');
imgv.node_data = vv.volt;

show_fem( imgv ); axis image
print_convert conformal1_02b.png '-density 125'
