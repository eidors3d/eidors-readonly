% Map elements onto mesh
% $Id: moving_ball01.m,v 1.1 2007-09-24 15:02:00 aadler Exp $

% Load models
load ng_mdl_16x1_coarse;
load ng_mdl_16x1_fine;
load ng_mdl_16x1_vfine;

[eptr, xyz]= mk_mesh_sample_array( ng_mdl_16x1_coarse, ...
                                   ng_mdl_16x1_fine, ...
                                   ng_mdl_16x1_vfine, ...
                                   128^3);
% This is very slow - save results
save ng_mdl_16x1_ptrs eptr xyz
