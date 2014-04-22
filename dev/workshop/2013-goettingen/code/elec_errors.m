clear
n_elecs = 32;
[stim,msel]= mk_stim_patterns(n_elecs,1,[0,5],[0,5],{'no_meas_current'}, 1);
fmdl = ng_mk_cyl_models([2 2 0.2] ,[n_elecs,1],[0.1]);
fmdl.stimulation =  stim;

opt.imgsz = [32 32];
opt.noise_figure = 0.5;
imdl = mk_GREIT_model(mk_image(fmdl,1), 0.25, [], opt);

vv= eidors_readdata('elec_errors.eit', 'LQ2');
vg= vv(msel,31:39); vh= mean(vg,2); % good
vg= vv(msel,61:69); vh= mean(vg,2); % bad

imgg = inv_solve(imdl, vh,vg);
subplot(211);show_slices(imgg);eidors_colourbar(imgg)

if 0
    imdl.meas_icov_rm_elecs.elec_list = [26,27];
    imdl.meas_icov_rm_elecs.exponent  = -1;
    imdl.meas_icov_rm_elecs.SNR       = 100;
    opt.noise_covar = meas_icov_rm_elecs(imdl);
else
   imdl.calc_reciproc_error.tau = -3e-3;
   opt.noise_covar = calc_reciproc_error( imdl, vg );
end
imdl = mk_GREIT_model(mk_image(fmdl,1), 0.25, [], opt);

%imgg.calc_colours.component = 'imag';
imgb = inv_solve(imdl, vh,vg);
subplot(212);show_slices(imgb); eidors_colourbar(imgb)