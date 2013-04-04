% Forward solvers $Id$

% Calculate a stimulation pattern
stim = mk_stim_patterns(19,1,[0,1],[0,1],{},1);

% Solve all voltage patterns
img_2.fwd_model.stimulation = stim;
img_2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img_2);

% Show first stim pattern
h1= subplot(221);
img_v = rmfield(img_2, 'elem_data');
img_v.node_data = vh.volt(:,1);
show_fem(img_v);

% Show 7th stim pattern
h2= subplot(222);
img_v = rmfield(img_2, 'elem_data');
img_v.node_data = vh.volt(:,7);
show_fem(img_v);

img_v.calc_colours.cb_shrink_move = [0.3,0.8,-0.02];
common_colourbar([h1,h2],img_v);
print_convert forward_solvers02a.png
