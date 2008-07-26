% run all tutorials in this directory
% $Id$

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
