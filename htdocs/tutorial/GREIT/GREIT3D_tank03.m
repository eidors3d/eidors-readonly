% STEP 3: Reconstruction model 
fmdl= ng_mk_cyl_models([4,1,.5],[16,1.5,2.5],[0.05]);
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(32,1,[0,5],[0,5],{},1);
fmdl = mdl_normalize(fmdl, 0);
[~,fmdl] = elec_rearrange([16,2],'square', fmdl);
