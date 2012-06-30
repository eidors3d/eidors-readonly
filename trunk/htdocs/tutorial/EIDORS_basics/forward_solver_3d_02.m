extra={'ball','solid ball = sphere(0,0.2,0.5;0.1);'};
fmdl= ng_mk_cyl_models([1,0.3,0.05],[nelec,ring_vert_pos],[0.1,0.05,0.02], extra);
fmdl.stimulation = stim;

img= mk_image(fmdl, conduct);
img.elem_data(fmdl.mat_idx{2}) = 0.1;

show_fem(img);
print_convert forward_solver_3d_02a.png '-density 75'
