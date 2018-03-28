% STEP 2: Simulate data
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(32,1,[0,5],[0,5],{},1);
img= mk_image(fmdl,1);
vh = fwd_solve(img);
for i=1:length(posns)-4;
   img= mk_image(fmdl,1);
   img.elem_data(fmdl.mat_idx{i+1}) = 2;
   vi{i} = fwd_solve(img);
end;
