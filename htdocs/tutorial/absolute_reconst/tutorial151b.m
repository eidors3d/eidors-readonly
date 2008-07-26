% Simulate data for difference and absolute imaging
% $Id$

% Create homogeneous simulation Model
backgnd= .02; 
sim_img.elem_data= backgnd*ones(size(simmdl.elems,1),1);
v_homg= fwd_solve(sim_img);

% Create target simulation Model
sim_img.elem_data(target1)= backgnd*2;
sim_img.elem_data(target2)= backgnd/4;
v_targ= fwd_solve(sim_img);

clf; subplot(211);
plot([v_homg.meas, v_targ.meas]);
print -r75 -dpng tutorial151b.png;
