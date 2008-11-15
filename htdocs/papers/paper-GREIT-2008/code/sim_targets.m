function sim_targets( savefile )
% simulate radial movement $Id$

params= [0.9,0.05,0.5,0.5]; %max_posn, targ_rad, z_0, z_t
stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}', {'no_meas_current'}, 1);

load ng_cyl_mdl.mat
fmdl.stimulation = stim_pat;

[vh,vi,xyzr_pt]= simulate_3d_movement(200, fmdl, params, @simulation_radmove);
%Change: mdl geometry at 90 deg; radius is 15
xyzr_pt= xyzr_pt([2,1,3,4],:)/15;

save(savefile, 'vh','vi','xyzr_pt');

return
% OLD CODE FOR 2D Simulations
   imdl = mk_common_model('f2d3c',16); fmdl= imdl.fwd_model;
   fmdl.stimulation = stim_pat;
   [vh,vi,xyzr_pt]= simulate_2d_movement(200, fmdl, params, @simulation_radmove);
   xyzr_pt = [-1,0,0;0,1,0;0,0,0;0,0,1]*xyzr_pt;
