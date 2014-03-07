function netgen_gen_models( numbers )

if nargin==0; numbers = [1:15,15.1,16,16.1,17:20]; end
for i=numbers(:)';
   do_sim(i);
end

function fmdl= do_sim( number )
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
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 img = mk_image( fmdl, 1);
 img.elem_data(fmdl.mat_idx{2}) = 1.1; 
 
 fmdl = img; % so that the code shows the image

  case 10;
 shape_str = ['solid top    = plane(0,0,0;0,0,1);\n' ...
              'solid mainobj= top and orthobrick(-3,-3,-2;3,3,0) -maxh=0.5;\n'];
 [elec_pos_x,elec_pos_y] = meshgrid(linspace( -1.5,1.5,5),linspace(-2,2,7));
 elec_pos = [  elec_pos_x(:), elec_pos_y(:), ones(size(elec_pos_x(:)))*[0,0,0,1] ];
 elec_shape=[0.2];
 elec_obj = 'top';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

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

   case 12;
shape_str = ['solid top     = ellipsoid(0,0,0; 0,0,1; 1,0,0; 0,1,0); \n' ...
    'solid mainobj= top and orthobrick(-2,-2,0;2,2,2) -maxh=0.1;\n'];
deg2rad = pi/180;
th = (-70:20:70)'*deg2rad;
 elec_pos = [0*th,sin(th),cos(th),0*th,sin(th),cos(th); ...
             sin(th),0*th,cos(th),sin(th),0*th,cos(th)];
 elec_shape=[0.05];
 elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 13;
 shape_str = ['solid cyl    = cylinder (0,0,0; 0,1,0; 1); \n', ...
              'solid bottom = plane(0, 0,0;0,-1,0);\n' ...
              'solid top    = plane(0,10,0;0, 1,0);\n' ...
              'solid cut1   = plane(0, 4,0;0,-1,0);\n' ...
              'solid cut2   = plane(0, 6,0;0, 1,0);\n' ...
              'solid ball   = cyl and cut1 and cut2;  tlo ball;\n' ...
              'solid mainobj= ( top and (not cut2) and cyl ) or ' ...
                      '(bottom      and (not cut1) and cyl ) -maxh=0.8;\n'];
 elec_pos = [ 0, 10,  0, 0,  1,  0; 
              0,  0,  0, 0, -1,  0]; 
 elec_shape=[1.0];
 elec_obj = {'top','bottom'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 fmdl = mk_image(fmdl,1); 
 fmdl.elem_data(fmdl.fwd_model.mat_idx{2}) = 1.1;

   case 14;
shape_str = ['solid cyl    = cylinder (0,0,0; 0,1,0; 1); \n', ...
             'solid bottom = plane(0, 0,0;0,-1,0);\n' ...
             'solid top    = plane(0,15,0;0, 1,0);\n' ...
             'solid cut1   = plane(0, 6,0;0,-1,0);\n' ...
             'solid cut2   = plane(0, 9,0;0, 1,0);\n' ...
             'solid rod   = cyl and (top and not cut2) or cyl and (bottom and not cut1);\n' ... 
             'tlo rod -maxh=0.8;\n',...
             'solid mainobj= orthobrick(-5,-5,-5;5,20,5) and not rod;\n'];

elec_pos = [0,15, 0, 0, 1, 0; %top end electrode
            0,13, 0,-1, 0, 0; %top cylinder electrode
            0, 0, 0, 0,-1, 0; %bot end electrode
            0, 2, 0, 1, 0, 0];%bot cylinder electrode
elec_shape=[6,0,0];
elec_obj = {'rod','rod','rod','rod'};
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
img  = mk_image(fmdl,1); 

% join electrodes
fmdl.electrode(3).nodes = unique(horzcat(fmdl.electrode(3:4).nodes));
fmdl.electrode(4) = [];
fmdl.electrode(1).nodes = unique(horzcat(fmdl.electrode(1:2).nodes));
fmdl.electrode(2) = [];

   case 15;
shape_str = [ ...
  'solid cyl     = cylinder (-2,0,0;-2,0,1; 0.1); \n', ...
  'solid borehole= cyl and plane (0,0,-2;0,0,-1); \n', ...
  'solid cel1= orthobrick(-10,-10,-1.1;10,10,-1.0) and cyl; tlo cel1 -maxh=0.5;\n' ...
  'solid cel2= orthobrick(-10,-10,-1.6;10,10,-1.5) and cyl; tlo cel2 -maxh=0.5;\n' ...
  'solid top    = plane(0,0,0;0,0,1);\n' ...
  'solid top_obj= top and orthobrick(-5,-5,-5;5,5,0) -maxh=0.5;\n' ...
  'solid mainobj= top_obj and not borehole;\n'];
[elec_pos_x,elec_pos_y] = meshgrid(linspace( -1.5,1.5,3),linspace(-2,2,4));
elec_pos = [  elec_pos_x(:), elec_pos_y(:), ones(size(elec_pos_x(:)))*[0,0,0,1] ];
elec_shape=[0.1];
elec_obj = repmat({'top'}, 1, size(elec_pos,1));
elec_pos(end+1,:) = [-2.0,0,-1.05,0,0,1]; elec_obj(end+1)   = {'cel1'};
elec_pos(end+1,:) = [-2.0,0,-1.55,0,0,1]; elec_obj(end+1)   = {'cel2'};
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 15.1;
% First run 15
crop_model([], inline('z<-2 | abs(x)>2.5 | abs(y)>2.5','x','y','z'));
axis([-2.5,2.5,-2.5,2.5,-2,0]);

   case 16;
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
 elec_shape=[0.4,1,0.05];
 elec_obj = {'ymin', 'ymin', 'ymin','ymax','xmin'};
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);

   case 16.1;
 fmdl = do_sim( 16 );
 fmdl = mdl2d_from3d(fmdl);

   case 17;
 shape_str = [ ...
  'solid cyl    = cylinder (0,0,0; 0,0,1; 1); \n', ...
  'solid mainobj= plane(0,0,0;0,0,-1)\n' ...
        'and  plane(0,0,2;0,0,1)\n' ...
        'and  cyl -maxh=0.3;\n' ...
  'solid in_elec = sphere(0,-1,1;0.2)' ...
        'and not    sphere(0,-1,1;0.15) -maxh=0.05;\n' ...
        'solid in_elec0= in_elec  and mainobj;\n' ...
        'tlo in_elec0 cyl;\n' ...
  'solid out_elec = sphere(0,-1,1;0.4)' ...
        'and not    sphere(0,-1,1;0.35) -maxh=0.05;\n' ...
        'solid out_elec0= out_elec  and mainobj;\n' ...
        'tlo out_elec0 cyl;\n'];
 % Find a background electrode (for all) first. This will stop errors in the next
 elec_pos = [  0, -1,   0, NaN,NaN,NaN;
               1,  0,   1,   1,  0,  0;
               0,  1, 1.2,   0,  1,  0;
               0, -1, 1.2, NaN,NaN,NaN;
               0, -1, 1.4, NaN,NaN,NaN];
 elec_shape=[0.1];
 elec_obj = 'cyl';
 fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
 % Throw away the first electrode (background). We don't need it
 fmdl.electrode = fmdl.electrode(2:end);

end

if ~exist('fmdl'); return; end

show_fem(fmdl);
if any(number==[11]); view(270,60); end
if any(number==[13]); view(-64,-13); end
if any(number==[14]); view(-111,21); end
if any(number==[15]); view(-10,20); end

if rem(number,1) ==0 
  prname = sprintf('netgen_gen_models%02d.png',floor(number));
else
  prname = sprintf('netgen_gen_models%02da.png',floor(number));
end
print_convert( prname, '-density 75');
