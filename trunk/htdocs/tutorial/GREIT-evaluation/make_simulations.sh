#!/bin/sh 
# Run matlab to make simulations

nohup \
matlab -nodisplay \
 -r "run ../../../eidors/startup;  tic; \
     params= [0.9,0.05,0.5,0.5]; \
     stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}', {'no_meas_current'}, 1); \
     load ng_tank_75; fmdl.stimulation = stim_pat; \
     [vh,vi,xyzr_pt]= simulate_3d_movement(200, fmdl, params, @simulation_radmove); \
     xyzr_pt= xyzr_pt([2,1,3,4],:)/15; \
     save sim_radmove_homog vh vi xyzr_pt; \
\
     toc; disp('FINISHED'); exit;"
