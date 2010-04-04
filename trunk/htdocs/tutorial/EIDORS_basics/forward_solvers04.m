% Forward solvers $Id$

% Calculate a stimulation pattern
stim = mk_stim_patterns(19,1,[0,1],[0,1],{},1);

% Solve all voltage patterns
img_1.fwd_model.stimulation = stim;
img_1.fwd_solve.get_all_meas = 1;
vh1= fwd_solve(img_1);

img_2.fwd_model.stimulation = stim;
img_2.fwd_solve.get_all_meas = 1;
vh2= fwd_solve(img_2);

img_v = rmfield(img_2, 'elem_data');

% Show homoeneous image
subplot(231);
img_v.node_data = vh1.volt(:,1);
show_fem(img_v);

% Show inhomoeneous image
subplot(232);
img_v.node_data = vh2.volt(:,1);
show_fem(img_v);

% Show difference image
subplot(233);
img_v.node_data = vh1.volt(:,1) - vh2.volt(:,1);
img_v.calc_colours.cb_shrink_move = [0.3,0.6,+0.03];
show_fem(img_v,1);

print_convert forward_solvers04a.png
