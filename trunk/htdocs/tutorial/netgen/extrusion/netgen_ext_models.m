function netgen_ext_models( numbers )

if nargin==0; numbers = 1:20; end
for i=numbers(:)';
   disp(i);
   do_sim(i);
end

function do_sim( number )
switch number


% Extruded 3D shape with 9 rectuangular electrodes
   case 1;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
fmdl = ng_mk_extruded_model({2,xy,1},[9,0,1],[0.05,0.4]);

% Extruded 3D shape with mesh refinement and 5 large circular electrodes
   case 2;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
fmdl = ng_mk_extruded_model({2,xy,1,0.1},[5,0,1,2],[0.10]);

% Extruded 3D shape mesh with boundary interpolation
   case 3;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
fmdl = ng_mk_extruded_model({1,xy,[4,50]},[8,0,0.3,0.6],[0.08]);

% 2D shape mesh boundary interpolation
   case 4;
   % had to change number of interpolation points from 50 to 47 and
   % electrodes to equidistant to avoid netgen playing up.
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
fmdl = ng_mk_extruded_model({0,xy,[4,47]},[16,1],[0.1]);

% 2D shape with object and specific electrode positions
   case 5;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.05;'};
elec_pos = [0, 0; 30,0;60,0;90,0];
fmdl = ng_mk_extruded_model({0,xy,[4,49]},elec_pos,[0.1],extra);
img = mk_image(fmdl,1);
img.elem_data( fmdl.mat_idx{2} ) = 1.1;
fmdl= img;

% 3D shape with two objects and no electrodes
   case 6;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
extra={'twoball','solid twoball = sphere(0.2,0.2,.7;0.2) or sphere(-0.2,-0.2,.4;0.3);'};
fmdl = ng_mk_extruded_model({1,xy,[4,20]},[0,0],[0.05],extra);
img = mk_image(fmdl,1);
img.elem_data( fmdl.mat_idx{2} ) = 1.1;
fmdl= img;

% 3D shape with custom electrodes
   case 7;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
elec_pos = [0, 0.1; 30,0.2; 60,0.3; 90,0.4;120,0.5];
extra={'cyl','solid cyl = cylinder(0.2,-0.2,0;0.2,-0.2,1;0.2) and orthobrick(-1,-1,0;1,1,1);'};
fmdl = ng_mk_extruded_model({1,xy,[4,47]},elec_pos,[0.05,0.2],extra);
img = mk_image(fmdl,1);
img.elem_data( fmdl.mat_idx{2} ) = 1.1;
fmdl= img;

% CAN'T MAKE IT WORK
   case 8.1;
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
xyi=[ -0.74 -0.74 -0.60 -0.36 -0.35 -0.55;
      -0.17 -0.51 -0.55 -0.44 -0.12  0.00]';
%[fmdl,mat_idx] = ng_mk_extruded_model({1,{xy,xyi},[4,20]},elec_pos,[0.05,0.2]);
fmdl = ng_mk_extruded_model({1,{xy,xyi},[0]},[0,0],[]);
img = mk_image(fmdl,1);
img.elem_data( fmdl.mat_idx{2} ) = 1.1;
fmdl= img;


   case 88;
xx=fliplr([
  -88.5777  -11.4887    4.6893   49.8963  122.7033  150.3033  195.5103 249.7573 ...
  258.8013  279.7393  304.9623  309.2443  322.0923  337.7963  340.6503 348.2633 ...
  357.3043  358.7333  361.5873  364.9183  365.3943  366.3453  366.3453 365.3943 ...
  362.5393  351.5943  343.5053  326.8513  299.2503  288.3073  264.9923 224.0703 ...
  206.4633  162.6833  106.5313   92.2543   57.5153    7.0733   -8.6297 -42.4167 ...
  -90.9547 -105.7057 -134.2577 -178.0367 -193.2647 -222.7687 -265.5957 -278.9197 ...
 -313.1817 -355.5337 -363.6237 -379.3267 -397.8857 -400.7407 -401.6927 -398.8377 ...
 -395.0307 -384.0867 -368.3837 -363.6247 -351.7277 -334.1217 -328.4117 -314.1357 ...
 -291.2947 -282.7297 -267.0257 -236.5707 -221.8187 -196.5977 -159.4807 -147.5837]);

yy=fliplr([
 -385.8513 -386.8033 -386.3273 -384.8993 -368.7193 -353.9673 -323.0363 -283.5403 ...
 -274.9743 -254.0363 -225.4843 -217.8703 -187.4153 -140.7813 -124.6013  -86.0573 ...
  -38.4703  -29.4273   -9.9173   21.0137   32.4347   53.3727   83.8257   93.3437 ...
  114.7587  149.0237  161.8717  187.5677  222.3037  231.3447  247.5237  267.5087 ...
  271.3177  277.0297  281.3127  279.4097  274.6507  273.2227  276.5547  284.6447 ...
  295.1127  297.4927  301.7757  304.1557  302.2537  297.4947  287.5017  282.2667 ...
  259.9017  225.6387  213.7427  185.6677  141.4127  125.2337   88.5917   34.8187 ...
   17.6897  -22.2803  -73.6723  -85.0923 -117.9263 -163.6083 -176.4573 -205.9613 ...
 -245.93 -256.4023 -275.4373 -304.9403 -315.4083 -332.0623 -352.0473 -355.3783]);
% Equidistant electrodes
fmdl = ng_mk_extruded_model({1,[xx',yy']/300,[4,40]},[16,1,0.5],[0.05]);

case 234918;
fmdl = ng_mk_extruded_model({2,{a,0.3*a},1},[16,0,1],[0.01]);
fmdl = ng_mk_extruded_model({2,{trunk/100, lung/100} ,[4,40]},[16,0,1],[0.1]);

% Thorax with two layers of electrodes
   case 8;
load CT2.mat;
fmdl = ng_mk_extruded_model({2,trunk/100 ,[4,50]},[6,0,0.5,1.5],[0.1,0.2]);


% Thorax with lung with 16 circular electrodes
   case 9;
load CT2.mat; lung = flipud(lung(1:3:end,:));
fmdl = ng_mk_extruded_model({1,{trunk/100, lung/100} ,[4,50]},[16,0,0.5],[0.1]);
img = mk_image(fmdl,1);
img.elem_data( fmdl.mat_idx{2} ) = 0.9;
fmdl= img;

end
if ~exist('fmdl'); return; end

show_fem(fmdl); view(170,20);
print_convert( ...
   sprintf('netgen_ext_models%02d_1.png',number), '-density 75');
view(0,90);
print_convert( ...
   sprintf('netgen_ext_models%02d_2.png',number), '-density 75');
