% Create simulation data $Id$

% 3D Model
imdl_3d= mk_common_model('n3r2',[16,2]);
show_fem(imdl_3d.fwd_model);

sim_img= mk_image( imdl_3d.fwd_model, 1);

% set homogeneous conductivity and simulate
homg_data=fwd_solve( sim_img );

% set inhomogeneous conductivity and simulate
sim_img.elem_data([390,391,393,396,402,478,479,480,484,486, ...
                   664,665,666,667,668,670,671,672,676,677, ...
                   678,755,760,761])= 1.15;
sim_img.elem_data([318,319,321,324,330,439,440,441,445,447, ...
                   592,593,594,595,596,598,599,600,604,605, ...
                   606,716,721,722])= 0.8;
inh_data=fwd_solve( sim_img );


slice_posn = [inf,inf,2.2,1,1; ...
              inf,inf,1.5,1,2;
              inf,inf,0.8,1,3];
show_slices(sim_img,slice_posn);

print -r75 -dpng basic_iterative01a.png
