% Map elements onto mesh
% $Id: moving_ball_16x2_01.m,v 1.1 2007-09-25 11:21:44 aadler Exp $

% Load models
load ng_mdl_16x2_coarse;
load ng_mdl_16x2_fine;
load ng_mdl_16x2_vfine;

[eptr, xyz]= mk_mesh_sample_array( ng_mdl_16x2_coarse, ...
                                   ng_mdl_16x2_fine, ...
                                   ng_mdl_16x2_vfine, ...
                                   128^3);
% This is very slow - save results
save ng_mdl_16x2_ptrs eptr xyz
