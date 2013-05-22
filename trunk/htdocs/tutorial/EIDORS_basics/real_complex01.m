imdl = mk_common_model('h2d1c',16);
img = mk_image(imdl,1);
vh = fwd_solve(img);

DeltaC1 = -0.1;    %Target 1 is non_conductive
DeltaC2 = 0+0.1i;  %Target 2 has + permittivity
target= mk_c2f_circ_mapping(img.fwd_model, [[0.5;0.0;0.1],[-0.5;0;0.1]]);
img.elem_data = 1+ DeltaC1*target(:,1) + DeltaC2*target(:,2) ;
vi = fwd_solve(img);
vi = add_noise(5,vi,vh);

img.calc_colours.component = 'real';
subplot(221); show_fem(img);
title 'real conductivity change'
print_convert real_complex01a.png 

img.calc_colours.component = 'imag';
subplot(221); show_fem(img);
title 'imag conductivity change'
print_convert real_complex01b.png 

