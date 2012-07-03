% Create pipe model $Id$

n_elec = 12;
stim= mk_stim_patterns( n_elec, 1, [0,3],[0,1],{},1);

fmdl= ng_mk_cyl_models(4,[n_elec,2],[0.2,0.5,0.04]);
fmdl.stimulation = stim;

clf; subplot(121);
show_fem(fmdl)
print_convert('pipe01a.png','-density 100')

show_fem(fmdl); view([0,0]);
print_convert('pipe01b.png','-density 90')
