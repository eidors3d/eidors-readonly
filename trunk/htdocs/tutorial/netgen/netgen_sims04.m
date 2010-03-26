% Netgen simulation $Id$

% Calculate stimulation pattern offset-4
img.fwd_model.stimulation = mk_stim_patterns(9,1,[0,4],[0,1],{},1);

% Set homogeneous phantom
img.elem_data(:) = 1;
vh = fwd_solve(img);
imgn.node_data = vh.volt(:,1);

show_fem(imgn); axis equal; print -dpng -r100 netgen_sims04a.png

% Set inhomogeneity very insulating
img.elem_data = 1 - 0.99*(ctr < 0.2^2);
vi = fwd_solve(img);
imgn.node_data = vi.volt(:,1);

show_fem(imgn); axis equal; print -dpng -r100 netgen_sims04b.png

% Difference in voltages
imgn.node_data = vh.volt(:,1) - vi.volt(:,1);

show_fem(imgn); axis equal; print -dpng -r100 netgen_sims04c.png
