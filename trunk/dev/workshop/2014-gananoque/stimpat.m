load montreal_data_1995
imdl = mk_common_model('b2c2',16);
fmdl = imdl.fwd_model;
show_fem(fmdl, [0, 1])
%% simulate data
[sp1, meas_sel] = mk_stim_patterns(16, 1,...
    [0 1], [0 1], {'rotate_meas'});
fmdl.stimulation = sp1;
himg = mk_image(fmdl,1);
vv = fwd_solve(himg);
clf
plot(vv.meas);

[sp2, meas_sel] = mk_stim_patterns(16, 1,...
    [0 1], [0 1], {'no_rotate_meas'});
fmdl.stimulation = sp2;
himg = mk_image(fmdl,1);
vv = fwd_solve(himg);
hold all
plot(vv.meas);

real_data = double(zc_h_demo3);
real_data = real_data(meas_sel);
plot(real_data/1e7)
axis tight