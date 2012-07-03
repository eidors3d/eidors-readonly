h= 10/2; w = 28; % Simulate tank of 30cm hight, 28cm width
stim= mk_stim_patterns(16,1,[0,1],[0,1],{},1); % Sheffield pattern
elec_sz = [0.2,0,0.05];   % electrode radius 0.5cm; 1cm diameter.

maxsz = .15;
fmdl = ng_mk_cyl_models([2*h,w/2,maxsz],[16,h],elec_sz); 
fmdl.stimulation = stim;

vr = fwd_solve(mk_image(fmdl,1)); % solve homogeneous model
vr.n_ne = [size(fmdl.nodes,1), size(fmdl.elems,1)];
