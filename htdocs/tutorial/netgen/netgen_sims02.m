% Netgen simulation $Id$

% Calculate stimulation pattern adjacent
stim = mk_stim_patterns(9,1,[0,1],[0,1],{},1);
img.fwd_model.stimulation = stim;

% Get all voltages so we can plot it
img.fwd_solve.get_all_meas = 1;

% Homogeneous model
img.elem_data(:) = 1;

vh = fwd_solve(img);

% Show Voltage for stim pattern #1
imgn = rmfield(img,'elem_data');
imgn.node_data = vh.volt(:,1);

show_fem(imgn);
print -dpng -r100 netgen_sims02a.png

% Show Voltage for stim pattern #2
imgn = rmfield(img,'elem_data');
imgn.node_data = vh.volt(:,2);

imgn.calc_colours.cb_shrink_move = [0.5,0.8,.02];

show_fem(imgn,1);
print -dpng -r100 netgen_sims02b.png
