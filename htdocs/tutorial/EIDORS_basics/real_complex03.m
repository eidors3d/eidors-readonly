% load data with complex measurements
load iirc_data_2006
vi= v_rotate(:,9); vh= v_reference;

imgr = inv_solve(imdl, vh, vi);

subplot(221); show_fem(imgr);
title 'real conductivity change'
print_convert real_complex03a.png 

imgr.calc_colours.component = 'imag';
subplot(221); show_fem(imgr);
title 'imag conductivity change'
print_convert real_complex03b.png 
