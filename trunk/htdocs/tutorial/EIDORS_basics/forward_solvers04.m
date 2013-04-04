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
h1= subplot(231);
img_v.node_data = vh1.volt(:,1);
show_fem(img_v);

% Show inhomoeneous image
h2= subplot(232);
img_v.node_data = vh2.volt(:,1);
show_fem(img_v);

% Show difference image
h3= subplot(233);
img_v.node_data = vh1.volt(:,1) - vh2.volt(:,1);
show_fem(img_v);

img_v.calc_colours.cb_shrink_move = [0.3,0.8,-0.05];
common_colourbar([h1,h2,h3],img_v);

print_convert forward_solvers04a.png
