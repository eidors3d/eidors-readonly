function img= inv_solve_dual_mesh( inv_model, data1)
% INV_SOLVE_DUAL_MESH using a coarse and fine mesh
% img= inv_solve_dual_mesh( inv_model, data1)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data object
%

% (C) 2005 David Stephenson. Licenced under the GPL Version 2
% $Id: inv_solve_dual_mesh.m,v 1.1 2005-12-06 19:40:11 aadler Exp $

voltage = data1.meas;

M_dense= inv_model.fwd_model;

[index_simp]=edge_refined_elem_mapper( mdl_coarse, mdl_dense);
c2d_elem= index_simp(:,1);
[index_vtx]=edge_refined_node_mapper( mdl_coarse, mdl_dense);
d2c_node= index_vtx(:,1);


% Load parameters
M_coarse= inv_model.inv_solve_dual_mesh.coarse_mdl;
mapper_func= inv_model.inv_solve_dual_mesh.mapper_func;

% Initial conductivity estimate
mat_ref_coarse = ones(size(M_coarse.elems,1),1);
mat_ref_coarse(:)= inv_model.jacobian_bkgnd.value;
mat_ref_dense= mat_ref_coarse( coarse2dense );

index_simp=edge_refined_elem_mapper(M_coarse,M_dense);

% create dense model
img_dense= eidors_obj('image', 'dense image');
img_dense.elem_data= mat_ref_dense;
img_dense.fwd_model= M_dense;
pp_dense= np_fwd_parameters( M_dense );

% create coarse model
img_coarse= eidors_obj('image', 'coarse image');
img_coarse.elem_data= mat_ref_coarse;
img_coarse.fwd_model= M_coarse;
pp_coarse= np_fwd_parameters( M_coarse );

%Set the tolerance for the forward solver
tol = 1e-5;

for iter= 1:inv_model.parameters.max_iterations;                     

   s_mat= calc_system_mat( M_dense, img_dense );
   V_dense = forward_solver( ...
              pp_dense.vtx, s_mat.E, pp_dense.I, tol, s_mat.perm);
   v_f_dense = m_3d_fields( ...
              pp_dense.vtx, pp_dense.n_elec, pp_dense.indH, ...
              s_mat.E, tol, pp_dense.gnd_ind);

    %% DENSE MESH TO COARSE MESH HERE for scalar potentials

   [V_coarse] = mapvd(V_dense,index_vtx, ...
                     length(M_dense.stimulation), ...
                     M_coarse.nodes, M_dense.nodes );

   [v_f_coarse] = mapvm(v_f_dense,index_vtx, ...
                     length(voltage), ...
                     M_coarse.nodes, M_dense.nodes );

   J_coarse= calc_jacobian( M_coarse, 



% create a data structure to return
img.name= 'solved by inv_solve_trunc_iterative';
img.elem_data = sol;
img.fwd_model= M_dense;

function voltH = elec_volts( Vfwd, fwd_model)
    p = np_fwd_parameters( fwd_model);
    Velec=Vfwd( p.n_node+(1:p.n_elec),:);
    voltH = zeros( p.n_meas, 1 );
    idx=0;
    for i=1:p.n_stim
       meas_pat= fwd_model.stimulation(i).meas_pattern;
       n_meas  = size(meas_pat,1);
       voltH( idx+(1:n_meas) ) = meas_pat*Velec(:,i);
       idx= idx+ n_meas;
    end

function [V_coarse] = mapvd(V_dense,index_vtx,n_stim_patterns,vtx_coarse,vtx_dense);

% This function maps the scaler potential field
% from the dense mesh to the coarse
% mesh using verticies extraction.

for j=1:n_stim_patterns

    for i=1:size(vtx_coarse,1);

        V_coarse(i,j)=V_dense((index_vtx(i,1)),j);
        
    end
    
end

V_coarse_elec=V_dense((size(vtx_dense,1)+1):(size(V_dense,1)),:);

V_coarse=[V_coarse;V_coarse_elec];

% This function maps the scaler potential field
% from the dense mesh to the coarse
% mesh using verticies extraction.
function [v_f_coarse] = mapvm(v_f_dense,index_vtx,n_measurements,vtx_coarse,vtx_dense);


for j=1:size(indH_dense,1);

    for i=1:size(vtx_coarse,1);

        v_f_coarse(i,j)=v_f_dense((index_vtx(i,1)),j);
        
    end
    
end

v_f_coarse_elec=v_f_dense((size(vtx_dense,1)+1):(size(v_f_dense,1)),:);

v_f_coarse=[v_f_coarse;v_f_coarse_elec];
