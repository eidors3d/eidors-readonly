xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';

%  Fourier desriptor
pp = fourier_fit(xy,length(xy)-1); % avoid overfitting

%  smoothened contour
ls = linspace(0,1,46); ls(end) = [];
xy = fourier_fit(pp,ls);

%  counter-clockwise contour
xy = flipud(xy);

fmdl = ng_mk_2d_model({xy,0.1},16,[0.1, 15]);
show_fem(fmdl,[0 1 0])