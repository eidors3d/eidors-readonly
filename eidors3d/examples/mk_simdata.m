% How to make simulation data using EIDORS3D
% $Id: mk_simdata.m,v 1.5 2005-06-04 19:00:11 aadler Exp $

% 
% Example 1: Create simple 16 electrode 2D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
n_elec= 16;
n_rings= 1;
params= mk_circ_tank(8, [], n_elec, n_rings); 
mdl_2d = eidors_obj('fwd_model', params);
mdl_2d.stimulation= mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}', ...
                            {'no_meas_current'}, 10);
mdl_2d.solve= 'aa_fwd_solve';
mdl_2d.system_mat= 'aa_calc_system_mat';
show_fem( mdl_2d );

% create homogeneous image + simulate data
mat= ones( size(mdl_2d.elems,1) ,1);

homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', mdl_2d );

homg_data=fwd_solve( mdl_2d , homg_img);

% 
% Example 1: Create simple 16 electrode 3D model
% 
% get parameters for model from mk_circ_tank
% param= mk_circ_tank(rings, levels, n_elec, n_planes )
params= mk_circ_tank(8, [-.5,0,.5], n_elec, n_rings); 
mdl_3d = eidors_obj('fwd_model', params);
