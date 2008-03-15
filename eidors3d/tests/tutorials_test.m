% Rerun all the tutorials to see if they work
% $Id: tutorials_test.m,v 1.2 2008-03-15 21:23:19 aadler Exp $

 cd ../../eidors-htdocs/tutorial

% EIDORS_basics
cd EIDORS_basics
clear; clf; close all; 
disp('tutorial010');
tutorial010a
tutorial010b
tutorial010c

clear; clf; close all; 
disp('tutorial110');
tutorial110a

clear; clf; close all; 
disp('tutorial120');
tutorial120a
tutorial120b

clear; clf; close all; 
disp('tutorial130');
tutorial130a
tutorial130b

clear; clf; close all; 
disp('backproj_solve');
backproj_solve01
backproj_solve02
clear; clf; close all; 
tutorial120a
backproj_solve03

% Find and trim all output images
!find -name '*.png' -exec convert -trim '{}' '{}' ';'

% adv_image_reconst
cd ../adv_image_reconst

disp('basic_iterative');
clear; clf; close all; 
basic_iterative01
basic_iterative02

disp('total_variation');
clear; clf; close all; 
total_variation01
total_variation02
total_variation03
total_variation04

disp('move2d');
clear; clf; close all; 
move_2d01
move_2d02

disp('move3d');
clear; clf; close all; 
move_3d01
move_3d02

disp('temporal_solver');
clear; clf; close all; 
temporal_solver01
temporal_solver02
temporal_solver03
temporal_solver04

disp('TV_hyperparams');
clear; clf; close all; 
TV_hyperparams01
TV_hyperparams02
name_base= 'tv_hp_00_img-a1=%3.1f-a2=%3.1f.jpg';
TV_hyperparams03
fid= fopen('TV-params-NSR=0.html','w');
TV_hyperparams04

% add noise
TV_hyperparams05
name_base= 'tv_hp_20_img-a1=%3.1f-a2=%3.1f.jpg';
TV_hyperparams03
fid= fopen('TV-params-NSR=.01.html','w');
TV_hyperparams04


% Don't know what to do - asked Camille
%vitro_2d01.m

% Find and trim all output images
!find -name '*.png' -exec convert -trim '{}' '{}' ';'
return


