% $Id$

img = mk_image(fmdl,1);
vh = fwd_solve(img);

select_fcn = inline('(x-1).^2+(y-1).^2+(z-1).^2<.03','x','y','z');
img.elem_data = 1 + 0.1*elem_select(img.fwd_model, select_fcn);
vi = fwd_solve(img);

show_fem(img); view(0,70);
print_convert mk_GREIT_matrix03a.png '-density 60'
