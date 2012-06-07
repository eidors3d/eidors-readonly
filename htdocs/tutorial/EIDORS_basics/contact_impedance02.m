fmdl.stimulation(1).stim_pattern = [0;1;0;-1];
fmdl.stimulation(1).meas_pattern = [0;1;0;-1]';
fmdl.solve =      @fwd_solve_1st_order;
fmdl.system_mat = @system_mat_1st_order;
fmdl.electrode(1).z_contact = 0.01;

img = mk_image(fmdl,1);
img.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img);
imgv= rmfield(img,'elem_data');
imgv.node_data = vh.volt;
show_fem(imgv);
axis([-1.1,1.1,-0.5,0.5]);

print_convert contact_impedance02a.png
