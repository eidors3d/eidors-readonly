% Create 3D model of a tunnel $Id$ 

% Simulation protocol. TODO: we need a geophysical stim protocol
stim = mk_stim_patterns(N_elec, 1, [0,4], [0,4], {'no_meas_current'},1);
fmdl.stimulation = stim;
cond_mdl = .1; % in S/m units
img = mk_image( fmdl, cond_mdl); 
vs_h = fwd_solve( img);

img.elem_data = cond_mdl*(1 + mk_c2f_circ_mapping(fmdl, [0;1.4;0;0.2]) );
vs_i = fwd_solve( img);

show_fem(img); ylim(2*[-1,1]); zlim(2*[-1,1]);

view(90,0); print_convert tunnelsim02a.png
