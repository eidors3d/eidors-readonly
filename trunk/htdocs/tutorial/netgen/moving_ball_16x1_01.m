% Map elements onto mesh
% $Id$

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
