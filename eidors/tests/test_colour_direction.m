% Test colour mapping
% $Id$

% 2D Test
imdl= mk_common_model('b2c',16);
s_img= eidors_obj('image','');
s_img.fwd_model= imdl.fwd_model;
s_img.elem_data= ones(size(s_img.fwd_model.elems,1),1);
v_homg = fwd_solve(s_img);
s_img.elem_data(1)= 1.1;
v_targ = fwd_solve(s_img);
imgr1= inv_solve(imdl,v_homg,v_targ);
slice1= calc_slices(imgr1);
subplot(121)
show_fem(imgr1);
disp('image should be red');

% 3D Test
imdl= mk_common_model('n3r2',[16 2]);
s_img= eidors_obj('image','');
s_img.fwd_model= imdl.fwd_model;
s_img.elem_data= ones(size(s_img.fwd_model.elems,1),1);
v_homg = fwd_solve(s_img);
s_img.elem_data(1)= 1.1;
v_targ = fwd_solve(s_img);
imgr2= inv_solve(imdl,v_homg,v_targ);
subplot(122)
show_fem(imgr2);
disp('image should be red');
