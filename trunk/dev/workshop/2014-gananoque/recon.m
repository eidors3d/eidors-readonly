vv = eidors_readdata('P04P-1016.get');

fmdl = ng_mk_cyl_models([0 1 .2], 16,[.1 0 .05]);
[fmdl.stimulation fmdl.meas_select] = ...
    mk_stim_patterns(16,1,'{ad}','{ad}');
%normalize
fmdl = mdl_normalize(fmdl,1);

imdl = select_imdl(fmdl,{'Basic GN dif'});%, 'Choose NF=1'});
imdl.hyperparameter.value = 3e-2;
imdl.inv_solve_diff_pdipm.norm_data  =2;
imdl.inv_solve_diff_pdipm.norm_image =1 ;
imdl = rmfield(imdl, 'RtR_prior');
imdl.R_prior = 'prior_TV';
imdl.solve = 'inv_solve_TV_pdipm';
imdl.inv_solve.calc_solution_error = 1;
rimg = inv_solve(imdl, vv(:,35), vv(:,45));
show_slices(rimg);
% plot(sum(rimg.elem_data,1))

