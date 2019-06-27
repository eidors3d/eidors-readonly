function netgen_cyl_models( numbers )

if nargin==0; numbers = 1:20; end
for i=numbers(:)';
   do_sim(i);
end
!for i in netgen_cyl_models*.png ; do convert -trim -depth 8 $i $i ; done

function do_sim( number )
switch number
% Simple 3D cylinder. Radius = 1. No electrodes
case 01;
    fmdl= ng_mk_cyl_models(3,[0],[]); 
show_fem(fmdl);

% Simple 2D cylinder. Radius = 2. Set minsize to refine
case 02;
    fmdl= ng_mk_cyl_models([0,2,.2],[0],[]); 
show_fem(fmdl);

% 3D cylinder. Radius = 1. 2 planes of 8 elecs with radius 0.1
case 03;
    fmdl= ng_mk_cyl_models(3,[8,1,2],[0.1]); 
show_fem(fmdl);

% 3D cylinder. Radius = 1. 6 circ elecs with elec refinement
case 04;
    fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0,0.05]); 
show_fem(fmdl);

% 3D cylinder. Radius = 1. 5 rect elecs with no refinement
case 05;
    fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0.3]); 
show_fem(fmdl);

% 2D cylinder. Radius = 1. 11 rect elecs with refinement
case 06;
    fmdl= ng_mk_cyl_models(0,[11],[0.2,0,0.05]); 
show_fem(fmdl);

% 2D cylinder. Radius = 1.5. Refined(0.1). 11 elecs with refinement
case 07;
    fmdl= ng_mk_cyl_models([0,1,0.1],[11],[0.2,0,0.02]); 
show_fem(fmdl);

% 2D cylinder. elecs at 0, 90 and 120 degrees
case 08;
    fmdl= ng_mk_cyl_models(0,[0;90;120],[0.2,0,0.03]); 
show_fem(fmdl);

% 2D cylinder. elecs at 0 (large,refined) and 90 (small) degrees
case 09;
    fmdl= ng_mk_cyl_models(0,[0;90],[0.4,0,0.01;0.1,0,0.1]); 
show_fem(fmdl);

% 3D cylinder. elecs at 0, 30, 60, 90 in planes
case 10;
    fmdl= ng_mk_cyl_models(3,[0,0.5;30,1;60,1.5;90,2.0],[0.2,0,0.1]); 
show_fem(fmdl);

% 3D cylinder. Various elecs at 0, 30, 60, 90 in planes
case 11;
    el_pos = [0,0.5;30,1;60,1.5;90,2.0];
    el_sz  = [0.2,0,0.1;0.1,0,0.05;0.2,0.2,0.02;0.2,0.4,0.5];
    fmdl= ng_mk_cyl_models(3,el_pos,el_sz); 
show_fem(fmdl);

% Simple 3D cylinder with a ball
case 12;
    extra={'ball','solid ball = sphere(0.5,0.5,2;0.4);'};
    fmdl= ng_mk_cyl_models(3,[0],[],extra); 
    img= mk_image(fmdl,1); img.elem_data(fmdl.mat_idx{2}) = 2;
show_fem(img);

% 3D cylinder with 8 electrodes and cube
case 13;
    extra={'cube','solid cube = orthobrick(0.5,0.5,0.5;0,0,1.5);'};
    fmdl= ng_mk_cyl_models(2,[8,0.5,1.5],[0.1],extra); 
    img= mk_image(fmdl,1); img.elem_data(fmdl.mat_idx{2}) = 2;
show_fem(img);

% 3D cylinder with inner cylinder
case 14;
    extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,1;1,1,2) -maxh=0.05;'};
    fmdl= ng_mk_cyl_models(3,[0],[],extra); 
show_fem(fmdl);

% 2D cylinder with 8 electrodes and hole
case 15;
    extra={'ball','solid ball = sphere(0.2,0.2,0;0.2) -maxh=0.05;'};
    fmdl= ng_mk_cyl_models(0,[8],[0.1,0,0.05],extra); 
show_fem(fmdl);

% 2D cylinder with 9 electrodes and inner cylinder
case 16;
    extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.03;'};
    fmdl= ng_mk_cyl_models(0,[9],[0.2,0,0.05],extra); 
    ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.2).^2 + (ctr(:,2)-0.2).^2;
    img= mk_image(fmdl, 1 + 0.1*(ctr<0.2^2));
show_fem(img);

% 2D cylinder with 16 electrodes and inner cylinder and box
case 17;
extra={'ballandbox', ...
       ['solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and ' ...
                     'orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;' ...
        'solid box = orthobrick(-0.3,-0.3,0;-0.1,0.1,0.05) -maxh=0.1;' ...
        'solid ballandbox = ball or box;']};

fmdl= ng_mk_cyl_models(0,[16],[0.1,0,0.05],extra); 
img = mk_image( fmdl, 1);
img.elem_data( fmdl.mat_idx{2} ) = 1.1;
show_fem(img);

% 3D cylinder with 8 electrodes and cube
case 18;
    ob= 'orthobrick(-1,-1,0;1,1,1.5)';
    el= 'ellipsoid(0,0,1.5; 0,0,1.25; 0.8,0,0; 0,0.8,0)';
%   sp= 'sphere(0, 0, 1.5; 0.8)';
    extra={'obj',['solid obj = ',ob,' and ',el,';']};
    fmdl= ng_mk_cyl_models(1.5,[8,0.5,1.0],[0.1],extra); 
    img= mk_image(fmdl,1); img.elem_data(fmdl.mat_idx{2}) = 2;
show_fem(img);

case 19;
    t1= 'solid obj1 = torus( 0, 0,  .8; 1, 0, 0; .4 ; .12);';
    t2= 'solid obj2 = torus( 0, 0, 1.2; 0, 1, 0; .4 ; .12);';
    extra={'obj1','obj2',[t1,t2]};
    fmdl= ng_mk_cyl_models(2,[6,.5],[0.3,0.2],extra);
    img= mk_image(fmdl,1);
    img.elem_data(fmdl.mat_idx{2}) = 1.5;
    img.elem_data(fmdl.mat_idx{3}) = 0.5;
show_fem(img);

case 20;
   extra={'ball_inside','ball_surface', [ ...
          'solid ball_inside  = sphere(-0.4,0,0.5;0.05);' ...
          'solid ball_surface = sphere( 0.4,0,1.0;0.05);' ...
          ]};
   fmdl= ng_mk_cyl_models(1,[8,.5],[.1],extra); 
   fmdl = mat_idx_to_electrode(fmdl, {2,3});
   show_fem(fmdl); view(3,12);
end

%print('-dpng','-r75', ...
%   sprintf('netgen_cyl_models%02d.png',number));
print_convert(  ...
    sprintf('netgen_cyl_models%02d.png',number), '-density 75');
