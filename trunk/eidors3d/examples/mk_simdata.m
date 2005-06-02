% How to make simulation data using EIDORS3D
% $Id: mk_simdata.m,v 1.1 2005-06-02 16:34:40 aadler Exp $

% get parameters for model fro mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
params= mk_circ_tank(8, [], 16, 1); 
sim_mdl = eidors_obj('fwd_model', params);
show_fem( sim_mdl );
