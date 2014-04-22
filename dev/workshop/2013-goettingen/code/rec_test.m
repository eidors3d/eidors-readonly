load montreal_data_1995;
vh = zc_h_demo4; vi=zc_demo4(:,1:4:16);
vh = zc_resp(:,34); vi = zc_resp(:,1:4:16);

fmdl = ng_mk_cyl_models(0,16,0.05);
[stim,mpat] = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
fmdl.stimulation = stim;
fmdl.meas_select = mpat;
%vh= vh(mpat); vi= vi(mpat,:);

imdl = select_imdl(fmdl, {'Basic GN dif'});
imdl = select_imdl(imdl, {'Choose NF=2'});
img = inv_solve(imdl, vh, vi);
subplot(411)
img.show_slices.img_cols = 4; show_slices(img);

imdl.RtR_prior = @prior_noser;
imdl = select_imdl(imdl, {'Choose NF=2'});
subplot(412)
img = inv_solve(imdl, vh, vi);
img.show_slices.img_cols = 4; show_slices(img);

imdl.RtR_prior = @prior_tikhonov;
imdl = select_imdl(imdl, {'Choose NF=2'});
subplot(413)
img = inv_solve(imdl, vh, vi);
img.show_slices.img_cols = 4; show_slices(img);

imdl = select_imdl(fmdl, {'TV solve dif'});
imdl.hyperparameter.value = 3;
subplot(414)
img = inv_solve(imdl, vh, vi);
img.show_slices.img_cols = 4; show_slices(img);




return





fmdl = ng_mk_cyl_models(0,16,0.05);

[fmdl.stimulation, fmdl.meas_select] = ...
    mk_stim_patterns(16,1,[0,2],[0,2],{},1);

subplot(411)
imdl = select_imdl(fmdl, {'Basic GN dif'});
img = inv_solve(imdl, vh, vi);
img.show_slices.img_cols = 4; show_slices(img);

[fmdl.stimulation, fmdl.meas_select] = ...
    mk_stim_patterns(16,1,[0,1],[0,1],{},1);

subplot(412)
imdl = select_imdl(fmdl, {'Basic GN dif'});
img = inv_solve(imdl, vh, vi);
img.show_slices.img_cols = 4; show_slices(img);


return
imdl.RtR_prior = @prior_noser;
imdl = select_imdl(imdl, {'Choose NF=2'});

img = inv_solve(imdl,zc_h_demo4, zc_demo4);


