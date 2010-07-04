% run all tutorials in this directory
% $Id$

clear; clf; close all; 
disp('centre_slice');
!7za x ~/docs/eidors/htdocs/data_contrib/netgen_moving_ball/ng_mdl_16x2_vfine.7z
centre_slice01
centre_slice02
centre_slice03
centre_slice04
centre_slice05
centre_slice06

clear; clf; close all; 
disp('dual_3d3d');
dual_3d3d

clear; clf; close all; 
disp('dual_model');
dual_model01
dual_model02
dual_model03
dual_model04
dual_model05

clear; clf; close all; 
disp('dual_partial');
dual_partial2d01
dual_partial2d02
dual_partial2d03
dual_partial2d04
dual_partial2d05

clear; clf; close all; 
disp('nodal_jacobian');
nodal_jacobian01
nodal_jacobian02
nodal_jacobian03

clear; clf; close all; 
disp('pipe');
pipe01
pipe02
pipe03
pipe04

clear; clf; close all; 
disp('square_mesh');
!7za x ../../data_contrib/netgen_moving_ball/ng_mdl_16x1_coarse.7z
!7za x ../../data_contrib/netgen_moving_ball/ng_mdl_16x1_vfine.7z
square_mesh01
square_mesh02
square_mesh03
square_mesh04
square_mesh05
square_mesh06

clear; clf; close all; 
disp 'nodal solver'
nodal_solve01
nodal_solve02
nodal_solve03

clear; clf; close all; 
disp 'rpi data'
rpi_data01
rpi_data02a
rpi_data02
rpi_data03
rpi_data04
rpi_data05
rpi_data06

clear; clf; close all; 
disp 'common models'
common_models01
common_models02
common_models03
common_models04

disp('two_and_half');
two_and_half_d01
two_and_half_d02
two_and_half_d03
