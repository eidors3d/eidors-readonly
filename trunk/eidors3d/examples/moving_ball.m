inv3d= eidors_obj('inv_model', 'EIT inverse');
inv3d.reconst_type= 'difference';
inv3d.jacobian_bkgnd.value = 1;
inv3d.fwd_model.misc.perm_sym= '{y}';
inv3d.hyperparameter.value = 1e-4;
inv3d.solve=            'np_inv_solve';
inv3d.R_prior.func=     'np_calc_image_prior';
inv3d.np_calc_image_prior.parameters= [3 1]; %  deg=1, w=1
inv3d.fwd_model= tank_mdls(1);

imgr= inv_solve(inv3d,tank_data(2),tank_data(1));

