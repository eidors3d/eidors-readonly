function netgen_gen_models( numbers )

if nargin==0; numbers = 1:20; end
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

   case 9;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid ball   = sphere(-1.25,0,-1;0.5); tlo ball;\n' ...
              'solid mainobj= top and orthobrick(-2,-1,-2;2,1,0) and not ball -maxh=0.5;\n'];
 elec_pos = linspace( -1.5,1.5,5)';
 elec_pos = [  elec_pos, elec_pos*[0,0,0,0], elec_pos*0+1];
 elec_shape=[0.3];
 elec_obj = 'top';
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 img = mk_image( fmdl, 1);
 img.elem_data(mat_idx{2}) = 1.1; 
 
 fmdl = img; % so that the code shows the image

  case 10;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-3,-3,-2;3,3,0) -maxh=0.5;\n'];
 [elec_pos_x,elec_pos_y] = meshgrid(linspace( -1.5,1.5,5),linspace(-2,2,7));
 elec_pos = [  elec_pos_x(:), elec_pos_y(:), ones(size(elec_pos_x(:)))*[0,0,0,1] ];
 elec_shape=[0.2];
 elec_obj = 'top';
 [fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 11;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,0,1; 1.0); \n', ...
              'solid tank   = orthobrick(-2,-2,0;2,2,0.4) and cyl; \n', ...
              'solid fish   = ellipsoid(0.2,0.2,0.2;0.2,0,0;0,0.1,0;0,0,0.1); tlo fish;\n', ...
              'solid mainobj= tank and not fish -maxh=0.3;\n'];
 n_elec = 7;
 th = linspace(0,2*pi,n_elec+1)'; th(end) = [];
 cs = [cos(th), sin(th)];
 elec_pos = [  cs, 0.2+0*th, cs, 0*th];
 elec_shape=[0.05];
 for i=1:n_elec; elec_obj{i} = 'cyl'; end
 i=i+1;elec_pos(i,:) = [ 0  ,0.2,0.2,-1,0,0]; elec_obj{i} = 'fish';
 i=i+1;elec_pos(i,:) = [ 0.4,0.2,0.2, 1,0,0]; elec_obj{i} = 'fish';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

end

if ~exist('fmdl'); return; end

show_fem(fmdl);
if any(number==[11]); view(270,60); end

print_convert( ...
   sprintf('netgen_gen_models%02d.png',number), '-density 75');
