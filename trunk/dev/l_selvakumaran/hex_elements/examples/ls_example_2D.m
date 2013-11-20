[fwd_mdl, mat_indices]= ls_gmsh_mk_fwd_model('mesh/trial_2D.msh', 'square', 'elec',[], 0.01);
fwd_mdl.system_mat='ls_system_mat_1st_order';
[stim, meas_sel]= mk_stim_patterns(4,1, '{ad}', '{ad}', [], 1);
fwd_mdl.stimulation= stim;
fwd_mdl.meas_select= meas_sel;
img=mk_image(fwd_mdl,1);

img.fwd_solve.get_all_meas = 1;
data = fwd_solve(img);


h1= subplot(121);
img_v = rmfield(img, 'elem_data');
img_v.node_data = data.volt(:,1);
show_fem(img_v);
h2= subplot(122);
img_v = rmfield(img, 'elem_data');
img_v.node_data = data.volt(:,2);
show_fem(img_v);
img_v.calc_colours.cb_shrink_move = [0.3,0.8,-0.02];
common_colourbar([h1,h2],img_v);

min(data.volt(:,2))
max(data.volt(:,2))

