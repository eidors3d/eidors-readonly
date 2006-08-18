% Explore Stimulation Patterns - and reconstruct images
% $Id: tutorial030c.m,v 1.1 2006-08-18 19:07:32 aadler Exp $

% Define an image
sim_img= eidors_obj('image', 'stimulation image');
sim_img.elem_data= ones( size(zigzag_mdl.elems,1) ,1);

% Simulate homogeneous measurements
sim_img.fwd_model= zigzag_mdl;
zigzag_data_h= fwd_solve( sim_img );

sim_img.fwd_model= planar_mdl;
planar_data_h= fwd_solve( sim_img );

sim_img.fwd_model= pl_ops_mdl;
pl_ops_data_h= fwd_solve( sim_img );

% Create targets in image
sim_img.elem_data([390,391,393,396,402,478,479,480,484,486, ...
                   664,665,666,667,668,670,671,672,676,677, ...
                   678,755,760,761])= 1.15;
sim_img.elem_data([318,319,321,324,330,439,440,441,445,447, ...
                   592,593,594,595,596,598,599,600,604,605, ...
                   606,716,721,722])= 0.8;

% Simulate inhomogeneous measurements
sim_img.fwd_model= zigzag_mdl;
zigzag_data_i= fwd_solve( sim_img );

sim_img.fwd_model= planar_mdl;
planar_data_i= fwd_solve( sim_img );

sim_img.fwd_model= pl_ops_mdl;
pl_ops_data_i= fwd_solve( sim_img );


%Add 25dB SNR noise to data
noise= std(zigzag_data_i.meas - zigzag_data_h.meas) ...
       / 10^(25/20) * randn(size(zigzag_data_h.meas));

zigzag_data_i.meas= zigzag_data_i.meas + noise;
planar_data_i.meas= planar_data_i.meas + noise;
pl_ops_data_i.meas= pl_ops_data_i.meas + noise;


%Reconstruct and show images
slices= [inf,inf,2.0,1,1; inf,inf,1.0,1,2];

subplot(131)
imdl_3d.fwd_model= zigzag_mdl;
img= inv_solve(imdl_3d, zigzag_data_i, zigzag_data_h);
show_slices(img, slices);

subplot(132)
imdl_3d.fwd_model= planar_mdl;
img= inv_solve(imdl_3d, planar_data_i, planar_data_h);
show_slices(img, slices);

subplot(133)
imdl_3d.fwd_model= pl_ops_mdl;
img= inv_solve(imdl_3d, pl_ops_data_i, pl_ops_data_h);
show_slices(img, slices);

print -r75 -dpng tutorial030c.png;
