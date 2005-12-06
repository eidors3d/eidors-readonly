% Example to show reconstructions from
% the Process Tomography 
% group from U.Manchester, UK
%
% (C) 2005 by Stephen Murphy. Licensed under GPL version 2.
% $Id: manchester_tomography.m,v 1.5 2005-12-06 21:14:36 aadler Exp $
function manchester_tomography( M_coarse, M_dense)

% load dual_mesh
% DS_coarse= set_fwd_model(vtx_coarse,simp_coarse,[],elec_coarse,zc,gnd_ind_coarse,Ib_coarse, [], []);
% DS_dense= set_fwd_model(vtx_dense,simp_dense,[],elec_dense,zc,gnd_ind_dense,Ib_dense, [], []);

example2( M_coarse, M_dense)

function example1
load 2planes_2rods_Opp
fwd_m1=eidors_obj('fwd_model',fwd_mdl);
fwd_m1=eidors_obj('fwd_model',fwd_mdl); 
fwd_m1.misc.perm_sym = '{y}';
img=eidors_obj('image','homog');


img.fwd_model= fwd_m1;
img.elem_data= mat_ref;
vh= fwd_solve(img);
%display_meas(fwd_m1)

%img.elem_data= ... 
%   drawsphere(fwd_m1.nodes,fwd_m1.elems,mat_ref,[0,0,.05],.05,.05);
%vi= fwd_solve(img);

if 0
fwd_m2= fwd_m1;
renumber=[1:2:16;18:2:32];renumber=renumber(:);
fwd_m2.electrode= fwd_m1.electrode(renumber);
fwd_m2.stimulation= mk_stim_patterns(16,1,'{ad}','{ad}',[],1);
img2=eidors_obj('image','homog');
img2.elem_data= mat_ref;
img2.fwd_model= fwd_m2;
vh2= fwd_solve(img2);
display_meas(fwd_m2)
end

fwd_m1.jacobian= 'np_calc_jacobian';

imdl= eidors_obj('inv_model','NP mdl');
if 1
    imdl.solve= 'np_inv_solve';
    imdl.solve= 'inv_solve_trunc_iterative';
    imdl.R_prior.func= 'np_calc_image_prior';
    imdl.hyperparameter.value = 1e-2;
else
    imdl.solve= 'ab_tv_diff_solve';
    imdl.R_prior.func= 'ab_calc_tv_prior';
    imdl.parameters.max_iterations= 2;
    imdl.hyperparameter.value = [1e-2, 1e-5];
end
imdl.jacobian_bkgnd.value = .01;
imdl.np_calc_image_prior.parameters= [3 1];
imdl.reconst_type= 'difference';
imdl.fwd_model= fwd_m1;

meas.type='data';
reference.type='data';
imr= inv_solve(imdl, meas, reference);
%imr= inv_solve(imdl, vi, vh);


show_slices(imr, linspace(.01,.09,4)'*[inf,inf,1])


function example2( M_coarse, M_dense)

stim= mk_stim_patterns(16,1,'{ad}','{ad}',[],1);
M_coarse.stimulation= stim;
M_dense.stimulation= stim;
%display_meas( M_dense);

n= size(M_dense.elems,1);
img2=eidors_obj('image','homog');
img2.elem_data= ones(n,1);
img2.fwd_model= M_dense;
vh= fwd_solve(img2);

cc= center_of_simps( M_dense);
[jnk,mat]= elems_in_cylinder(cc,[0,.03,.02],.01,ones(n,1),1.2);
img2.elem_data= mat;

show_fem(img2);
vi= fwd_solve(img2);

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
