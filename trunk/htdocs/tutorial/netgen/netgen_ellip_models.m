function netgen_ellip_models( numbers )

if nargin==0; numbers = 1:9; end
for i=numbers(:)';
   do_sim(i);
end

function do_sim( number )
switch number
% Simple 3D ellipse. Major, minor axes = [1.5, 0.8]. No electrodes
case 01;
    fmdl= ng_mk_ellip_models([1,1.5,0.8],[0],[]); 
show_fem(fmdl);

% Simple 2D cylinder. Axes = [1.5,0.8]. Refined to 0.1
case 02;
    fmdl= ng_mk_ellip_models([0,1.5,0.8,0.1],[],[]); 
show_fem(fmdl);

% 3D cylinder. Axes = [1.5,0.8]. 2 planes of 8 elecs with radius 0.1
case 03;
    fmdl= ng_mk_ellip_models([1,1.2,0.8],[8,0.3,0.7],[0.1]); 
show_fem(fmdl);

% 3D cylinder. Axes= [1.3,1] = 1. 7 rect elecs with no refinement
case 04;
    fmdl= ng_mk_ellip_models([3,1.3],[7,1],[0.2,0.3]); 
show_fem(fmdl);

% 2D cylinder. Axes = [1.2,0.8]. 11 rect elecs with refinement. Rotated 45 degrees
case 05;
    fmdl= ng_mk_ellip_models([0,1.2,0.8],[11],[0.2,0,0.05]); 
th = 45* [2*pi/360];
fmdl.nodes = fmdl.nodes*[cos(th),sin(th);-sin(th),cos(th)];
show_fem(fmdl);


% 2D cylinder. elecs at 0, 90 and 120 degrees
case 06;
    fmdl= ng_mk_ellip_models([0,1.2,0.8],[0;90;120],[0.2,0,0.03]); 
show_fem(fmdl);


% 3D cylinder. Various elecs at 0, 30, 60, 90 in planes
case 07;
    el_pos = [0,0.5;30,1;60,1.5;90,2.0];
    el_sz  = [0.2,0,0.1;0.1,0,0.05;0.2,0.2,0.02;0.2,0.4,0.5];
    fmdl= ng_mk_ellip_models([3,0.8,1.2],el_pos,el_sz); 
show_fem(fmdl);

% Simple 3D cylinder with a ball
case 08;
    extra={'ball','solid ball = sphere(0.5,0.5,1;0.4);'};
    fmdl= ng_mk_ellip_models([2,1.2,0.8],[8,1],[0.1],extra); 
    img= mk_image(fmdl, 1);
    img.elem_data(fmdl.mat_idx{2}) = 2;
show_fem(img);

case 09;
    b1 = 'solid ball1= sphere(0.5, 0.5,1;0.2);';
    b2 = 'solid ball2= sphere(0.5,-0.5,1;0.2);';
    extra = {'ball1','ball2',[b1,b2]};
    [fmdl,mat_idx]= ng_mk_ellip_models([2,1.2,0.8],[8,1],[0.1],extra);
    img= mk_image(fmdl, 1);
    img.elem_data(mat_idx{2}) = 2; 
    img.elem_data(mat_idx{3}) = 0.5;
show_fem(img);


end


print_convert( ...
   sprintf('netgen_ellip_models%02d.png',number), '-density 75');
