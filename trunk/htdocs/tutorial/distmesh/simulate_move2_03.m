% simulate homogeneous $Id$

% stimulation pattern: adjacent
stim_pat= mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},1);

smdl.stimulation= stim_pat;
himg= mk_image(smdl, 1);

vh= fwd_solve(himg);
