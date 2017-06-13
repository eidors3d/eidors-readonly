%Alternative Gauss Newton solver, changing prior at each iteration
imdl.solve = @inv_solve_abs_GN_prior;
imdl.hyperparameter.value = .01;
img   = inv_solve(imdl, v_i);     img.calc_colours   = cc;
img_n = inv_solve(imdl, v_i_n);  img_n.calc_colours = cc;
vr_agn=fwd_solve(img); vr_agn_n=fwd_solve(img_n);

clf; show_3d_slices(img,0.6,0.3,0.3); eidors_colourbar(img);
print_convert compare_3D_abs_GN_03a.png '-density 75'

clf; show_3d_slices(img_n,0.6,0.3,0.3);
eidors_colourbar(setfield(img_n,'calc_colours',cb));
print_convert compare_3D_abs_GN_03b.png '-density 75'
