inv3d= mk_common_model('n3r2',16);
inv3d.jacobian_bkgnd.value = 1;
inv3d.hyperparameter.value = 1e-1;
inv3d.solve=            'np_inv_solve';
inv3d.R_prior.func=     'np_calc_image_prior';
inv3d.np_calc_image_prior.parameters= [3 1]; %  deg=1, w=1
inv3d.fwd_model.stimulation= mk_stim_patterns( 16,2,'{op}','{ad}',{'meas_current'});
inv3d.fwd_model.nodes(:,[1,2])= inv3d.fwd_model.nodes(:,[1,2])*1.5;

 imgr= inv_solve(inv3d,tank_data(30),tank_data(1));

 show_slices(imgr,linspace(.1,2.9,9)'*[inf,inf,1])
