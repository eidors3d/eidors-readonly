% Example to show reconstructions from
% the Process Tomography 
% group from U.Manchester, UK
%
% manchester_tomography( example_no)
%
% example_no == 1
%    show reconstructions from 2planes / 2rods
%    Reconstruction = inv_solve_trunc_iterative
%
% example_no == 2
%    show reconstructions from 2planes / 2rods
%    Reconstruction = np_inv_solve
%
% (C) 2005 by Stephen Murphy. Licensed under GPL version 2.
% $Id: manchester_tomography.m,v 1.7 2005-12-07 15:09:11 aadler Exp $
function manchester_tomography( example_no)

switch example_no
    case 1,
        example_diff_morozov_reconst;

    case 2,
        example_diff_np_reconst;

    case 3,
        example_diff_tv_reconst

    case 10,
        example2
example2( M_coarse, M_dense)
end

% output fwd_model, vinh, vhom
function [fwd_m1,meas_2rod,meas_ref]= twoplane_mdl;
    load 2planes_Op_drive;
    fwd_m1=eidors_obj('fwd_model',mdl_2plane);
    fwd_m1.misc.perm_sym = '{y}';
    fwd_m1.solve=    'np_fwd_solve';
    fwd_m1.jacobian= 'ms_calc_jacobian';
    fwd_m1.system_mat= 'np_calc_system_mat';
    meas_ref.type='data';
    meas_2rod.type='data';

function [DS_coarse, DS_dense]= coarse_dense_mdl
    load dual_mesh;
    DS_coarse= set_fwd_model(vtx_coarse,simp_coarse,[], ...
              elec_coarse,zc,gnd_ind_coarse,Ib_coarse, [], []);
    DS_dense= set_fwd_model(vtx_dense,simp_dense,[], ...
              elec_dense,zc,gnd_ind_dense,Ib_dense, [], []);
    stim= mk_stim_patterns(16,1,'{ad}','{ad}',[],1);
    M_coarse.stimulation= stim;
    M_dense.stimulation= stim;


function example_diff_morozov_reconst
    [fwd_m1,vi,vh]= twoplane_mdl;

    imdl= eidors_obj('inv_model','NP mdl');
    imdl.solve= 'np_inv_solve';
    imdl.solve= 'inv_solve_trunc_iterative';
    imdl.R_prior.func= 'np_calc_image_prior';
    imdl.jacobian_bkgnd.value = .01;
    imdl.np_calc_image_prior.parameters= [3 1];
    imdl.reconst_type= 'difference';
    imdl.fwd_model= fwd_m1;

    imr= inv_solve(imdl, vi, vh);
    show_slices(imr, linspace(.01,.09,4)'*[inf,inf,1])

function example_diff_np_reconst
    [fwd_m1,vi,vh]= twoplane_mdl;

    imdl= eidors_obj('inv_model','NP mdl');
    imdl.solve= 'np_inv_solve';
    imdl.R_prior.func= 'np_calc_image_prior';
    imdl.hyperparameter.value = 1e-4;

    imdl.jacobian_bkgnd.value = .01;
    imdl.np_calc_image_prior.parameters= [3 1];
    imdl.reconst_type= 'difference';
    imdl.fwd_model= fwd_m1;

    imr= inv_solve(imdl, vi, vh);
    show_slices(imr, linspace(.01,.09,4)'*[inf,inf,1])

function example_diff_tv_reconst
    [fwd_m1,vi,vh]= twoplane_mdl;

    imdl= eidors_obj('inv_model','NP mdl');
    imdl.solve= 'ab_tv_diff_solve';
    imdl.R_prior.func= 'ab_calc_tv_prior';
    imdl.parameters.max_iterations= 1;
    imdl.hyperparameter.value = [1e-2,1e-5];

    imdl.jacobian_bkgnd.value = .01;
    imdl.np_calc_image_prior.parameters= [3 1];
    imdl.reconst_type= 'difference';
    imdl.fwd_model= fwd_m1;

    imr= inv_solve(imdl, vi, vh);
    show_slices(imr, linspace(.01,.09,4)'*[inf,inf,1])


function example_diff_sim_reconst
    load 2planes_Op_drive;
    fwd_m1=eidors_obj('fwd_model',mdl_2plane);
    fwd_m1.misc.perm_sym = '{y}';
    fwd_m1.solve=    'np_fwd_solve';
    fwd_m1.jacobian= 'np_calc_jacobian';
    fwd_m1.system_mat= 'np_calc_system_mat';

function [vi,vh] = sim_inhomg( mdl);
    n= size(mdl.elems,1);
    img2=eidors_obj('image','homog');
    img2.elem_data= ones(n,1);
    img2.fwd_model= mdl;
    vh= fwd_solve(img2);

    cc= center_of_simps( mdl);
    [jnk,mat]= elems_in_cylinder(cc,[0,.03,.02],.01,ones(n,1),1.2);
    [jnk,mat]= elems_in_cylinder(cc,[0,-.03,-.02],.01,mat,0.8);
    img2.elem_data= mat;

%   show_fem(img2);
    vi= fwd_solve(img2);

function [vi,vh,mdl_z]= zig_zag_from_2rings(mdl_2r)
    mdl_z= mdl_2r;
    renumber=[1:2:16;18:2:32];
    renumber=renumber(:);
    mdl_z.electrode= mdl_2r.electrode(renumber);
    mdl_z.stimulation= mk_stim_patterns(16,1,'{ad}','{ad}',[],1);

    img_z=eidors_obj('image','homog');
    img_z.fwd_model= mdl_z;

    n= size(mdl.elems,1);
    mat= ones(n,1);
    img_z.elem_data= mat;
    vh= fwd_solve(img_z);

    cc= center_of_simps( mdl_z);
    [jnk,mat]= elems_in_cylinder(cc,[0,.03,.02],.01,mat,1.2);
    [jnk,mat]= elems_in_cylinder(cc,[0,-.03,-.02],.01,mat,0.8);
    img_z.elem_data= mat;
    vi= fwd_solve(img_z);
%   display_meas(fwd_m2)

function example2( M_coarse, M_dense)
    [M_coarse, M_dense]= coarse_dense_mdl;

imdl= mk_common_model('b2c',16);
imdl.hyperparameter.value= 1e-2;
imgr= inv_solve(imdl, vi,vh);
show_fem(imgr);


imdl2= eidors_obj('inv_model', 'DS dual mesh');
imdl2.solve= 'inv_solve_dual_mesh';
imdl2.hyperparameter.value = 1e-4;
imdl2.R_prior.func= 'np_calc_image_prior';
imdl2.np_calc_image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
imdl2.jacobian_bkgnd.value= 1;
imdl2.reconst_type= 'static';
imdl2.fwd_model= M_dense;
imdl2.parameters.max_iterations= 5;

imdl2.inv_solve_dual_mesh.coarse_mdl= M_coarse;
imdl2.inv_solve_dual_mesh.mapper_func= 'edge_refined_elem_mapper';

imgr= inv_solve(imdl2, vi);
keyboard
