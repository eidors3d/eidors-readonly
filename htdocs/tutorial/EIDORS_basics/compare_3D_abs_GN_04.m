%Constrained Gauss Newton solver
imdl.solve = @inv_solve_gn;
imdl.hyperparameter.value = hp;
% limit conductivity to be greater than 0 with a log parametrization
imdl.inv_solve_gn.elem_working = 'log_conductivity';

img   = inv_solve(imdl, v_i);    img.calc_colours   = cc;
img_n = inv_solve(imdl, v_i_n);  img_n.calc_colours = cc;

clf; show_3d_slices(img,  0.6,0.3,0.3); eidors_colourbar(img);
print_convert compare_3D_abs_GN_04a.png '-density 75'

clf; show_3d_slices(img_n,0.6,0.3,0.3); eidors_colourbar(img_n);
print_convert compare_3D_abs_GN_04b.png '-density 75'
