 shape_str = ['solid top    = plane( 0, 0, 0; 0, 0, 1);\n' ...
              'solid bot    = plane( 0, 0,-1; 0, 0,-1);\n' ...
              'solid xmax   = plane( 3, 0, 0; 1, 0, 0);\n' ...
              'solid xmin   = plane(-3, 0, 0;-1, 0, 0);\n' ...
              'solid ymax   = plane( 0, 2, 0; 0, 2, 0);\n' ...
              'solid ymin   = plane( 0,-2, 0; 0,-2, 0);\n' ...
              'solid mainobj= top and bot and xmax and xmin and ymax and ymin;'];
 elec_pos = [  1, -2,  0,   0,  1,  0;
               0, -2,  0,   0,  1,  0;
              -1, -2,  0,   0,  1,  0;
              -1,  2,  0,   0, -1,  0;
              -3,  1,  0,  -1,  0,  0];
 elec_shape=[0.1,2];
 elec_obj = {'ymin', 'ymin', 'ymin','ymax','xmin'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
subplot(211);
 show_fem( fmdl );
subplot(212);
 fmdl2 = mdl2d_from3d(fmdl);
 show_fem( fmdl2,[0,1]);
 return
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-2,-2,-2;2,2,0);\n'];
 elec_pos = [  1,  0,  0,   0,  0,  1;
               0,  0,  0,   0,  0,  1;
              -1,  0,  0,   0,  0,  1];
 elec_shape=[0.1];
 elec_obj = 'top';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 show_fem( fmdl );
 fmdl2 = mdl2d_from3d(fmdl);
 show_fem( fmdl2);
disp(1)

% Run current directories
