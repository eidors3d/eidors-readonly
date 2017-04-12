% test c2f
rimg = mk_image(cmdl,1);
rimg.fwd_model = imdl.fwd_model;
rimg.rec_model = imdl.rec_model;
rimg.elem_data = 1 + ...
   elem_select( rimg.rec_model, '(x-0.3).^2+(y+0.1).^2<0.15^2' ) - ...
   elem_select( rimg.rec_model, '(x+0.4).^2+(y-0.2).^2<0.2^2' )*0.5;
fimg = mk_image(fmdl,1);
fimg.elem_data = 1 + ...
   elem_select( rimg.fwd_model, '(x-0.3).^2+(y+0.1).^2<0.15^2' ) - ...
   elem_select( rimg.fwd_model, '(x+0.4).^2+(y-0.2).^2<0.2^2' )*0.5;
clf; % figure 1
subplot(131); show_fem(mk_image(imdl.rec_model,rimg.elem_data)); axis square; axis off;
title('coarse');
subplot(132); show_fem(mk_image(imdl.fwd_model,c2f*rimg.elem_data)); axis square; axis off;
title('coarse-to-fine');
subplot(133); show_fem(fimg); axis square; axis off; title('fwd fine');
print_convert two_and_half_d02a.png '-density 75';

% Simulate data - homogeneous
himg = mk_image(fmdl,1);
vh = fwd_solve_2p5d_1st_order(himg);
% Simulate data - inhomogeneous
vi = fwd_solve_2p5d_1st_order(fimg);
