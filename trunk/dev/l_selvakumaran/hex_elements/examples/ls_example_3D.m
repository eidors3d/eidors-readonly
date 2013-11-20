[fwd_mdl, mat_indices]= ls_gmsh_mk_fwd_model('mesh/trial_3D.msh', 'plate', 'elec',[], 0.01);
[stim, meas_sel]= mk_stim_patterns(4,1, '{ad}', '{ad}', [], 1);
fwd_mdl.stimulation= stim;
fwd_mdl.meas_select= meas_sel;
img=mk_image(fwd_mdl,1);

img.fwd_solve.get_all_meas = 1;
data = fwd_solve(img);

min(data.volt(:,2))
max(data.volt(:,2))

