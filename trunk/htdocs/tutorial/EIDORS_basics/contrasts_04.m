sz = 5;
for ellipse_x = [0.5,1,2];

   shape_str = [ ...
    sprintf('solid ec = ellipticcylinder(0,0,0;%f,0,0;0,%f,0);\n', ...
        0.5*ellipse_x, 0.5/ellipse_x ) ...
    sprintf('solid left  = plane(-%f,0,0;-1,0,0);\n',sz) ...
    sprintf('solid right = plane( %f,0,0; 1,0,0);\n',sz) ...
    sprintf('solid brick = orthobrick(-%f,-2,-0.2;%f,2,0);\n',sz+1,sz+1) ...
    'solid cyl = ec and brick; tlo cyl;\n' ...
    'solid mainobj= left and right and (not cyl) and brick;\n'];
   elec_pos = [ sz,  0,  0,   0,  0,  1;
               -sz,  0,  0,   0,  0, -1];
   elec_shape=[2.0];
   elec_obj = {'left','right'};
   fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
   fmdl = mdl2d_from3d(fmdl);

   fmdl.stimulation = stim_meas_list([1,2,1,2]);
   img = mk_image(fmdl,1);
   img.fwd_model.solve =      @fwd_solve_1st_order;
   img.fwd_model.system_mat = @system_mat_1st_order;
   [img.fwd_model.electrode(:).z_contact] = deal(1000); % Large
   img.fwd_solve.get_all_meas = 1;
   for contrast = linspace( -2,2,5);
      img.elem_data( fmdl.mat_idx{1} ) = exp(contrast);
      vh = fwd_solve(img);
      contrasts_03;
   end
end

