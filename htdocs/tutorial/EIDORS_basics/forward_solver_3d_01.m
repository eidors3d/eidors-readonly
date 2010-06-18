nelec= 8; nrings= 2;
ring_vert_pos = [0.2, 0.5]; 
fmdl= ng_mk_cyl_models([1,0.3,0.05],[nelec,ring_vert_pos],[0.1,0.05,0.02]);

stim = mk_stim_patterns(nelec,nrings,[0,1],[0,1],{'meas_current'},1);
fmdl.stimulation = stim;

conduct = 1;
img = mk_image( fmdl, conduct ); 

show_fem(img);
print_convert forward_solver_3d_01a.png '-density 75'
