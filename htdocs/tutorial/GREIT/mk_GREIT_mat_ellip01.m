% $Id$
n_elecs = 16;
stim =  mk_stim_patterns(n_elecs,1,[0,1],[0,1],{'no_meas_current'}, 1);

extra={'lungs','solid lungs = sphere(0.9,0.1,1;0.6) or sphere(-0.9,0.1,1;0.6);'};
[fmdl,midx] = ng_mk_ellip_models([2, 2,1.4,0.2] ,[n_elecs,1],[0.1], extra);
fmdl.stimulation =  stim;

img = mk_image(fmdl,1); % Homogeneous background
vh = fwd_solve(img);
img.elem_data(midx{2}) = 0.5; % Lung regions
vi = fwd_solve(img);

show_fem(img); view(0,70);
print_convert mk_GREIT_mat_ellip01a.png '-density 60'
