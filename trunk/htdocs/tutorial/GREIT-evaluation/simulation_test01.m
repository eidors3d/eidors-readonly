load ng_mdl_16x1_fine;
params= [0.9,0.05,0.5,0.5]; %max_posn, targ_rad, z_0, z_t
ng_mdl_16x1_fine.stimulation = mk_stim_patterns(16,1,'{ad}','{ad}', ...
           {'no_meas_current'}, 1);
[vh,vi,xyzr_pt]= simulate_3d_movement(200, ...
        ng_mdl_16x1_fine, params, @simulation_radmove);
save sim_radmove_homog vh vi xyzr_pt
