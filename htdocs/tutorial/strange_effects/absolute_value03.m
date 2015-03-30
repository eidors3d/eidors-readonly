% create inverse model
imdl = mk_common_model('c2c2',16);
imdl.fwd_model.stimulation = stim;
imdl.fwd_model.meas_select = msel;
imdl.hyperparameter.value = 0.10;

% reconstruct normally
show_fem( inv_solve( imdl, vh, vi));
print_convert absolute_value03a.png

% reconstruct abs
vha = abs(vh);
via = abs(vi);
show_fem( inv_solve( imdl, vha, via));
print_convert absolute_value03b.png
