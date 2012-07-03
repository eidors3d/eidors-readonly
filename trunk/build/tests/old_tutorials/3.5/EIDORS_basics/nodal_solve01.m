% Simulate Target $Id$

imdl= mk_common_model('d2d1c',16);
img= mk_image(imdl);

vh = fwd_solve(img); %Homogeneous
select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.1^2','x','y','z');
img.elem_data = 1 + elem_select(img.fwd_model, select_fcn);
vi = fwd_solve(img); %Homogeneous

subplot(221);
show_fem(img);
print_convert nodal_solve01a.png

plot([vh.meas, vi.meas]); axis tight
print_convert nodal_solve01b.png
