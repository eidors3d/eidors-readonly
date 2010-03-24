% Simple model with three electrodes $Id$
elec_pts{1} = [1,0]; % Point electrode
elec_pts{2} = [0,1;sin(0.2),cos(0.2)]; % Complete electrode between points
elec_pts{3} = [0.5,0.5]; % Point internal electode

subplot(221);
fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.05] );
show_fem(fmdl)
hold on; plot(0.5,0.5,'o','Color',[0,0.5,0]); hold off

subplot(222);
fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.02] );
show_fem(fmdl)
hold on; plot(0.5,0.5,'o','Color',[0,0.5,0]); hold off

print -dpng -r125 circ_2d_model01.png
