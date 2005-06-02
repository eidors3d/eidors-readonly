% How to make simulation data using EIDORS3D
% $Id: mk_simdata.m,v 1.2 2005-06-02 17:13:49 aadler Exp $

% 
% Example 1: Create simple 16 electrode 2D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
params= mk_circ_tank(8, [], 16, 1); 
sim_mdl = eidors_obj('fwd_model', params);
show_fem( sim_mdl );
