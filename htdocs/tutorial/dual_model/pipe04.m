cmdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( fmdl, cmdl);

% modify
imdl.fwd_model.coarse2fine = c2f;
imdl.hyperparameter.value = .01;

imgc= inv_solve(imdl, vh, vi);
imgc.calc_colours.ref_level= 0;

show_fem(imgc);
print_convert pipe04a.png '-density 75';
