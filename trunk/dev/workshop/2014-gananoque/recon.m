vv = eidors_readdata('P04P-1016.get');

fmdl = ng_mk_cyl_models([0 1 .1], 16,[.1 0 .01]);
[fmdl.stimulation fmdl.meas_select] = ...
    mk_stim_patterns(16,1,'{ad}','{ad}');
%normalize
fmdl = mdl_normalize(fmdl,1);

imdl = select_imdl(fmdl,{'Basic GN dif'});
imdl.hyperparameter.value = 1e-4;
imdl.inv_solve.calc_solution_error = 0;
rimg = inv_solve(imdl, vv(:,35), vv(:,45));
show_slices(rimg);
% plot(sum(rimg.elem_data,1))

