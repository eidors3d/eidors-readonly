switch 3.13;
   case 1;
imdl= mk_common_model('a2c2',8); fmdl = imdl.fwd_model;
   case 1.1;
imdl= mk_common_model('a3cr',8); fmdl = imdl.fwd_model;
   case 2;
imdl= mk_common_model('f2c2',8); fmdl = imdl.fwd_model;
   case 2.1;
imdl= mk_common_model('h2c2',8); fmdl = imdl.fwd_model;
   case 3;
fmdl= ng_mk_cyl_models(1,[7,.3,.7],[0.1]); 
fmdl.stimulation = mk_stim_patterns(14,1,[0,1],[0,1],{},1);
   case 3.1;
fmdl= ng_mk_cyl_models(1,[13,.3],[0.03]); 
fmdl.stimulation = mk_stim_patterns(13,1,[0,1],[0,1],{},1);
   case 3.11;
fmdl= ng_mk_cyl_models(0,[11,.3],[0.03]); 
fmdl.stimulation = mk_stim_patterns(11,1,[0,1],[0,1],{},1);
   case 3.12;
fmdl= ng_mk_cyl_models([0,1,.01],[15,.3],[0.05]); 
fmdl.stimulation = mk_stim_patterns(15,1,[0,1],[0,1],{},1);
   case 3.13;
fmdl= ng_mk_cyl_models([.5,1,.03],[7,.3],[0.05]); 
fmdl.stimulation = mk_stim_patterns(7,1,[0,1],[0,1],{},1);
   case 3.2;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-3,-3,-2;3,3,0) -maxh=0.3;\n'];
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
disp([size(J0) size(J1) norm(J0(:) - J1(:))]);

