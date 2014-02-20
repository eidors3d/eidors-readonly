if 1
imdl= mk_common_model('a2c2',8); fmdl = imdl.fwd_model;
else
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-3,-3,-2;3,3,0) -maxh=0.1;\n'];
 [elec_pos_x,elec_pos_y] = meshgrid(linspace( -1.5,1.5,5),linspace(-2,2,7));
 elec_pos = [  elec_pos_x(:), elec_pos_y(:), ones(size(elec_pos_x(:)))*[0,0,0,1] ];
 elec_shape=[0.2];
 elec_obj = 'top';
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
end
img= mk_image(fmdl,1);
show_fem(img);
J0= jacobian_adjoint(img.fwd_model,img);
J1= test_calc_jacobian(img.fwd_model,img);
disp(norm(J0(:) - J1(:)));
