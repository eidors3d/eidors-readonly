% Forward solvers $Id$


% Solve all voltage patterns
img_2.fwd_solve.get_all_meas = 1;
img_2.fwd_model.mdl_slice_mapper.npx = 64;
img_2.fwd_model.mdl_slice_mapper.npy = 64;

% Show [0-3] stim pattern
subplot(221);
stim = mk_stim_patterns(19,1,[0,3],[0,1],{},1);
img_2.fwd_model.stimulation = stim;
vh = fwd_solve(img_2);
show_current(img_2,vh.volt(:,1));
axis([-1,1,-1,1]);

% Show [2-9] stim pattern
subplot(222);
stim = mk_stim_patterns(19,1,[0,7],[0,1],{},1);
img_2.fwd_model.stimulation = stim;
vh = fwd_solve(img_2);
show_current(img_2,vh.volt(:,3));
axis([-1,1,-1,1]);

print -dpng -r125 forward_solvers05a.png
