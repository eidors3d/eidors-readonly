% Create fwd model
el_pos = [190,0.5;170,0.5];
extra = {'cube',['solid cube = orthobrick(-0.2 ,-0.97,0.6;0.2,0,0.7) or ' ...
                              'orthobrick( 0.15,-0.97,0.4;0.2,0,0.7) ;']};
fmdl= ng_mk_cyl_models([1,1,.05],el_pos,[0.05,0,0.05],extra); 

% Solve fwd model
fmdl.stimulation(1).stim_pattern = [1;-1];
fmdl.stimulation(1).meas_pattern = [1,-1]; % dummy
img = mk_image(fmdl,1);
img.elem_data(fmdl.mat_idx{2}) = 10000;
img.fwd_solve.get_all_meas = 1;
vh=fwd_solve(img);

clf;show_fem(img);
print_convert view_3D_surf04a.jpg '-density 75'
view(0,0);
print_convert view_3D_surf04b.jpg '-density 75'

