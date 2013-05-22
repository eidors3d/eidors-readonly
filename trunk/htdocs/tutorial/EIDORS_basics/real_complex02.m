% Create reconstruction model + solve
imdl = mk_common_model('b2d1c',16);
imgr = inv_solve(imdl, vh, vi);

imgr.calc_colours.clim = 0.02;
subplot(221); show_fem(imgr);
title 'real conductivity change'
print_convert real_complex02a.png 

imgr.calc_colours.component = 'imag';
subplot(221); show_fem(imgr);
title 'imag conductivity change'
print_convert real_complex02b.png 
