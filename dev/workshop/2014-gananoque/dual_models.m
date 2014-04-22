fmdl = mk_library_model('neonate_16el_lungs');
[fmdl.stimulation fmdl.meas_select] = ...
    mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl = mdl_normalize(fmdl,1);
imdl = select_imdl(fmdl, {'Basic GN dif'});
imdl = mk_pixel_slice(imdl,[inf inf 0.5]);
%%
imdl.hyperparameter.value = 1e-4;
vv = eidors_readdata('P04P-1016.get');
rimg = inv_solve(imdl, vv(:,35), vv(:,45));
show_slices(rimg);