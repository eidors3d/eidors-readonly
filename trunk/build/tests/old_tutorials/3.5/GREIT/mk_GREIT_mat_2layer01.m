% $Id$
n_elecs = 8;
layers = [0.8,1.2];
stim =  mk_stim_patterns(n_elecs*length(layers),1,[0,1],[0,1], ...
             {'no_meas_current'}, 1);

extra={'lungs','solid lungs = sphere(0.9,0.1,1.65;0.3) or sphere(-0.9,0.1,1;0.3);'};
[fmdl,midx] = ng_mk_ellip_models([2, 2,1.4,0.2] ,[n_elecs,layers],[0.1], extra);
fmdl.stimulation =  stim;

img = mk_image(fmdl,1); % Homogeneous background
vh = fwd_solve(img);
img.elem_data(midx{2}) = 0.5; % Lung regions
vi = fwd_solve(img);

show_fem(img,[0,1]); view(0,70);
print_convert mk_GREIT_mat_2layer01a.png '-density 60'
view(0,10);
print_convert mk_GREIT_mat_2layer01b.png '-density 60'
