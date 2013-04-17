%Constrained Gauss Newton solver
imdl.solve = @inv_solve_abs_GN;
imdl.hyperparameter.value = hp;
imdl.inv_solve_abs_GN.min_value = 0.9;
imdl.inv_solve_abs_GN.max_value = 1.1;

img   = inv_solve(imdl, v_i);
img_n = inv_solve(imdl, v_i_n);

clf; show_3d_slices(img,0.6,0.3,0.3);
eidors_colourbar(setfield(img,'calc_colours',cb));
print_convert compare_3D_abs_GN_04a.png '-density 75'

clf; show_3d_slices(img_n,0.6,0.3,0.3);
eidors_colourbar(setfield(img_n,'calc_colours',cb));
print_convert compare_3D_abs_GN_04b.png '-density 75'
