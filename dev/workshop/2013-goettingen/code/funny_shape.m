function mdl = funny_shape
xy = ginput;
% save xy xy
% load xy
cla
plot(xy(:,1),xy(:,2));
hold on
pp = fourier_fit(xy,5);
frac = linspace(0,1,21); frac(end) = [];
xy = fourier_fit(pp,frac);
plot(xy(:,1),xy(:,2),'r');
hold off

mdl = ng_mk_2d_model({xy,0.1},-16);
subplot(212)
show_fem(mdl,[0 1 0])