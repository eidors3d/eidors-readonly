[stim,mpat]= mk_stim_patterns(16,1,[0,1],[0,1],{},1);
vh = zc_h_demo3(mpat); vi= zc_demo3(mpat,1);


[vh, vi, p] = face('large');
[vh, vi] = face('small');

vi = add_noise(5,vi,vh);

fmdl = funny_shape;
fmdl = ng_mk_cyl_models([0,1,.1],16,.05);
imdl1 = mk_common_model('c2c2',16);
fmdl1 = imdl1.fwd_model;

fmdl1.stimulation = stim;
fmdl.stimulation = stim;

imdl = select_imdl(fmdl, {'Basic GN dif'});
imdl1 = select_imdl(fmdl1, {'Basic GN dif'});
%imdl.RtR_prior = @prior_tikhonov;
% imdl.tutorial210_cheat_tikhonov.cheat_weight = .5;
% imdl.RtR_prior = @tutorial210_cheat_tikhonov;

imdl1.hyperparameter.value = .01;
imdl.hyperparameter.value = .01;
%imdl = select_imdl(imdl, {'Choose NF=2'});

subplot(211)
% imdl.tutorial210_cheat_tikhonov.cheat_elements = p.sad;
imgr = inv_solve(imdl1, vh, vi);
show_fem(imgr);

subplot(212)
% imdl.tutorial210_cheat_tikhonov.cheat_elements = p.happy;
imgr = inv_solve(imdl, vh, vi);
show_fem(imgr);
