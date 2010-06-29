function netgen_ellip_models( numbers )

if nargin==0; numbers = 1:8; end
for i=numbers(:)';
   do_sim(i);
end

function do_sim( number )
switch number
   case 1;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid p1= plane(0,0,0;0,0.5,-1);\n' ...
              'solid p2= plane(0,0,2;0.5,0, 1);\n' ...
              'solid mainobj= p1 and p2 and cyl -maxh=0.3;\n'];
 elec_pos = []; elec_shape = []; elec_obj = {};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 2;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= plane(0,0,0;0,0,-1)\n' ...
                        'and  plane(0,0,1.5;0.2,0.2,1)\n' ...
                        'and  cyl -maxh=0.3;\n'];
 elec_pos = [ -1,  0,  1,   1,  0,  0;
               0, -1,1.2,   0,  1,  0]; 
 elec_shape=[0.1,0,0;0.2,0.3,0.02];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 3;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= orthobrick(-2,-2,0;2,2,2) and cyl -maxh=0.3;\n'];
 th = linspace(0,2*pi,15)'; th(end) = [];
 cs = [cos(th), sin(th)];
 elec_pos = [  cs, th/2/pi + 0.5, cs, 0*th];
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 4;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid mainobj= orthobrick(-2,-2,0;2,2,2) and cyl -maxh=0.3;\n'];
 th = linspace(0,2*pi,15)'; th(end) = [];
 cs = [cos(th), sin(th)];
 elec_pos = [  cs, th/2/pi + 0.5, cs, 0*th];
 elec_shape=[0.1*th/2/pi + 0.05];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 5;
 shape_str = ['solid cyl    = cylinder (0,0,0; 1,0,0; 1); \n', ...
              'solid mainobj= plane(0,0,0;-1,0,0)\n' ...
                        'and  plane(2,0,0;1,0,0)\n' ...
                        'and  plane(0,0,-0.5;0,0,-1)\n' ...
                        'and  cyl -maxh=0.3;\n'];
 elec_pos = [  1,  0,  1,   0,  0,  1;
             1.2,  1,  0,   0,  1,  0]; 
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 6;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
              'solid bottom = plane(0,0,0;0,0,-1);\n' ...
              'solid top    = plane(0,0,2;0,0,1);\n' ...
              'solid mainobj= top and bottom and cyl -maxh=0.3;\n'];
 elec_pos = [  1,  0,  1,   1,  0,  0;
               0,  1,1.2,   0,  1,  0;
               0.8,  0,  0, 0,  0, -1]; 
 elec_shape=[0.1];
 elec_obj = {'cyl','cyl','bottom'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 7;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-2,-2,-2;2,2,0);\n'];
 elec_pos = [  1,  0,  0,   0,  0,  1;
               0,  0,  0,   0,  0,  1;
              -1,  0,  0,   0,  0,  1];
 elec_shape=[0.1];
 elec_obj = 'top';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 8;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid cyl    = ellipticcylinder(0,0,0;2.5,0,0;0,1,0);\n' ...
              'solid mainobj= top and cyl and orthobrick(-2,-2,-2;2,2,0);\n'];
 elec_pos = [  1,  0,  0,   0,  0,  1;
               0,  0,  0,   0,  0,  1;
              -1,  0,  0,   0,  0,  1;
               1, -1,-1.2,  0, -1,  0;
               0, -1,-1.0,  0, -1,  0;
              -1, -1,-0.8,  0, -1,  0];
 elec_shape=[0.1];
 elec_obj = {'top','top','top','cyl','cyl','cyl'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);


end

show_fem(fmdl);

print_convert( ...
   sprintf('netgen_gen_models%02d.png',number), '-density 75');
