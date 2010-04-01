% Forward solvers $Id$

% Calculate a stimulation pattern
stim = mk_stim_patterns(19,1,[0,9],[0,1],{},1);

% Solve all voltage patterns
img_2.fwd_model.stimulation = stim;
img_2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img_2);

% Show first stim pattern
subplot(221);
img_v = rmfield(img_2, 'elem_data');
img_v.node_data = vh.volt(:,1);
show_fem(img_v);

% Show 7th stim pattern
subplot(222);
img_v = rmfield(img_2, 'elem_data');
img_v.node_data = vh.volt(:,7);
img_v.calc_colours.cb_shrink_move = [0.3,0.6,+0.03];
show_fem(img_v,1);

print -dpng -r125 forward_solvers03a.png
