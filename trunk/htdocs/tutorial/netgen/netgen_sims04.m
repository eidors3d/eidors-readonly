% Netgen simulation $Id$

% Calculate stimulation pattern offset-4
img.fwd_model.stimulation = mk_stim_patterns(9,1,[0,4],[0,1],{},1);

% Set homogeneous phantom
img.elem_data(:) = 1;
vh = fwd_solve(img);
imgn.node_data = vh.volt(:,1);
imgn.calc_colours.cb_shrink_move = [0.3,0.7,0.02];
imgn.calc_colours.ref_level = 0;

clf; subplot(221);
show_fem(imgn,1); axis equal; axis tight;
print_convert netgen_sims04a.png '-density 100'

% Set inhomogeneity very insulating
img.elem_data = 1 - 0.99*(ctr < 0.2^2);
vi = fwd_solve(img);
imgn.node_data = vi.volt(:,1);

clf; subplot(221);
show_fem(imgn,1); axis equal; axis tight;
print_convert netgen_sims04b.png '-density 100'

% Difference in voltages
imgn.node_data = vh.volt(:,1) - vi.volt(:,1);

clf; subplot(221);
show_fem(imgn,1); axis equal; axis tight;
print_convert netgen_sims04c.png '-density 100'
