imdl = mk_common_model('b2c2',16);
fmdl = imdl.fwd_model;

SKIP = 5;

sp   = mk_stim_patterns(16,1,[0 SKIP],[0 1]);
fmdl.stimulation = sp;
himg = mk_image(fmdl,1);
vv   = fwd_solve(himg);
clf
plot(vv.meas(1:32));
