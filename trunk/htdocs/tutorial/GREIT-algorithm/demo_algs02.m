% Simulated objects
imb=  mk_common_model('c2c',16);

img= calc_jacobian_bkgnd( imb );
vv= fwd_solve( img ); v(2).vh= vv.meas;
img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])= 1.1;
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])= 1.1;
vv= fwd_solve( img ); v(2).vi= vv.meas;

subplot(221);show_fem(img);
axis square; axis off
print -r75 -dpng demo_algs02a.png;
