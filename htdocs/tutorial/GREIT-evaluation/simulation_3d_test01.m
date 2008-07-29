% simulate radial movement $Id$
use_3d_model = 0;

if exist('sim_radmove_homog.mat','file')
   load sim_radmove_homog.mat vh vi xyzr_pt
else
   params= [0.9,0.05,0.5,0.5]; %max_posn, targ_rad, z_0, z_t
   stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}', {'no_meas_current'}, 1);
   if use_3d_model;
%     load ng_mdl_16x1_fine; fmdl= ng_mdl_16x1_fine;
      load ng_tank; fmdl.stimulation = stim_pat;
      [vh,vi,xyzr_pt]= simulate_3d_movement(200, fmdl, params, @simulation_radmove);
      xyzr_pt= xyzr_pt([2,1,3,4],:)/15; %Change: mdl geometry at 90 deg; radius is 15
   else;
      imdl = mk_common_model('f2d3c',16); fmdl= imdl.fwd_model;
      fmdl.stimulation = stim_pat;
      [vh,vi,xyzr_pt]= simulate_2d_movement(200, fmdl, params, @simulation_radmove);
      xyzr_pt = [-1,0,0;0,1,0;0,0,0;0,0,1]*xyzr_pt;
   end

   save sim_radmove_homog vh vi xyzr_pt
end
