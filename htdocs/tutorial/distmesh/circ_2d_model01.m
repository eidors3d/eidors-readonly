% Simple model with three electrodes $Id$
elec_pts{1} = [1,0]; % Point electrode
elec_pts{2} = [0,1;sin(0.2),cos(0.2)]; % Complete electrode between points
elec_pts{3} = [0.5,0.5]; % Point internal electode

fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.05] );
subplot(221); show_fem(fmdl); axis equal

fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.02] );
subplot(222); show_fem(fmdl); axis equal

print_convert circ_2d_model01.png '-density 125'
