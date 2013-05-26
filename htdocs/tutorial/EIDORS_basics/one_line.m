% $Id$
load montreal_data_1995
show_slices(inv_solve(mk_common_model('d2c2',16),zc_h_demo4,zc_demo4))
print_convert one_line01a.png '-density 75'