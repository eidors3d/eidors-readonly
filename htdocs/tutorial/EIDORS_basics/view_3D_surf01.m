% Create fwd model
el_pos = [190,0.5;170,0.5];
fmdl= ng_mk_cyl_models(1,el_pos,[0.05,0,0.05]); 

% Solve fwd model
fmdl.stimulation(1).stim_pattern = [1;-1];
fmdl.stimulation(1).meas_pattern = [1,-1]; % dummy
img = mk_image(fmdl,1);
img.fwd_solve.get_all_meas = 1;
vh=fwd_solve(img);

clf;show_fem(fmdl);
print_convert view_3D_surf01a.png '-density 75'

