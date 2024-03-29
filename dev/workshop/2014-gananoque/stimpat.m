clf
data = eidors_readdata('P04P-1016.get');
hold all

imdl = mk_common_model('b2c2',16);
fmdl = imdl.fwd_model;
% show_fem(fmdl, [0, 1])
%% simulate data
[sp1, meas_sel] = mk_stim_patterns(16, 1,...
    [0 1], [0 1], {'rotate_meas'});
fmdl.stimulation = sp1;
himg = mk_image(fmdl,1);
vv = fwd_solve(himg);
plot(vv.meas);

[sp2, meas_sel] = mk_stim_patterns(16, 1,...
    [0 1], [0 1], {'no_rotate_meas'});
fmdl.stimulation = sp2;
himg = mk_image(fmdl,1);
vv = fwd_solve(himg);
plot(vv.meas);

real_data = data(meas_sel);
plot(real_data/1e3)
axis tight