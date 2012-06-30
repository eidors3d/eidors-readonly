% Set stimulation (elecs 1,2 are the fish)
fmdl.stimulation(1).stim_pattern = [1;-1;zeros(n_elec,1)];
fmdl.stimulation(1).meas_pattern = sparse([ ...
    0, 0, 1, 0, 0, 0,-1, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0,-1, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0,-1, 0;
    0, 0, 0, 0, 0, 1, 0, 0, 0,-1]);

img = mk_image( fmdl, 1);
img.elem_data( fmdl.mat_idx{1} ) = 1.1;  % fish is conductive
img.elem_data( fmdl.mat_idx{3} ) = 0.001;% target is non-conductive
