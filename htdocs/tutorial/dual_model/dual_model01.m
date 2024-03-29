% Simulate data $Id$
imdl= mk_common_model('c2c2',16);
img= mk_image(imdl);
vh=fwd_solve(img);
idx=[365,328,292,259,227,198,170,145,121];
img.elem_data([idx,idx+1,101,81])=1.1;
vi=fwd_solve(img);

show_fem(img);
print_convert('dual_model01a.png', '-density 60');
