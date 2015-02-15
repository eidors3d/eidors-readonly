img.fwd_model.solve =      @fwd_solve_1st_order;
img.fwd_model.system_mat = @system_mat_1st_order;
[img.fwd_model.electrode(:).z_contact] = deal(1000); % Large

img.fwd_solve.get_all_meas = 1;
vv = fwd_solve(img);
imgv= rmfield(img,'elem_data');
imgv.node_data = vv.volt;
imgv.calc_colours.ref_level = mean(vv.volt);
show_fem(imgv);

print_convert contrasts_02a.png
