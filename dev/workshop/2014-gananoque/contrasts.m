imdl = mk_common_model('b2c2',16);
fmdl = imdl.fwd_model;
img = mk_image(fmdl,1);
show_fem(img,[0 1 1]);
vh = fwd_solve(img);
select_fcn = inline('(x-0.2).^2+(y-0.5).^2<0.2^2','x','y','z');
es = elem_select(fmdl, select_fcn);
img.elem_data = 1 + es*0.5;
show_fem(img,[0 1 1]);
vi = fwd_solve(img);
plot(vh.meas - vi.meas);

rimg = inv_solve(imdl,vh,vi);
show_fem(rimg);