% Netgen simulation $Id$

% Calculate stimulation pattern adjacent
img.fwd_model.stimulation = mk_stim_patterns(9,1,[0,1],[0,1],{},1);
vh = fwd_solve(img);
imgn.node_data = vh.volt(:,1);

show_fem(imgn); axis equal; print_convert netgen_sims03a.png '-density 100'

% Calculate stimulation pattern offset-2
img.fwd_model.stimulation = mk_stim_patterns(9,1,[0,2],[0,1],{},1);
vh = fwd_solve(img);
imgn.node_data = vh.volt(:,1);

show_fem(imgn); axis equal; print_convert netgen_sims03b.png '-density 100'

% Calculate stimulation pattern offset-3
img.fwd_model.stimulation = mk_stim_patterns(9,1,[0,3],[0,1],{},1);
vh = fwd_solve(img);
imgn.node_data = vh.volt(:,1);

show_fem(imgn); axis equal; print_convert netgen_sims03c.png '-density 100'

% Calculate stimulation pattern offset-4
img.fwd_model.stimulation = mk_stim_patterns(9,1,[0,4],[0,1],{},1);
vh = fwd_solve(img);
imgn.node_data = vh.volt(:,1);

show_fem(imgn); axis equal; print_convert netgen_sims03d.png '-density 100'
