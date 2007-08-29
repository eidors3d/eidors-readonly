function img= inv_solve_dual_mesh( inv_model, voltage)
% INV_SOLVE_DUAL_MESH using a coarse and fine mesh
% img= inv_solve_dual_mesh( inv_model, data1)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data object
%

% (C) 2005 David Stephenson. License: GPL version 2 or version 3
% $Id: inv_solve_dual_mesh.m,v 1.5 2007-08-29 09:15:31 aadler Exp $

M_dense= inv_model.fwd_model;
% Load parameters
M_coarse= inv_model.inv_solve_dual_mesh.coarse_mdl;
%mapper_func= inv_model.inv_solve_dual_mesh.mapper_func;

[index_simp]=edge_refined_elem_mapper( M_coarse, M_dense);
c2d_elem= index_simp(:,1);
[index_vtx]=edge_refined_node_mapper( M_coarse, M_dense);
d2c_node= index_vtx(:,1);

%FIXME: create parameters
tol_for= 1e-5;
tol_inv= 1e-5;

tfac= calc_hyperparameter(inv_model);
IM_coarse= inv_model; IM_coarse.fwd_model= M_coarse;
RtR= calc_RtR_prior( IM_coarse );

% Initial conductivity estimate
mat_ref_coarse = ones(size(M_coarse.elems,1),1);
mat_ref_coarse(:)= inv_model.jacobian_bkgnd.value;
mat_ref_dense= mat_ref_coarse( c2d_elem );
% Initial estimate - homogeneous background coarse mesh
sol_upd_dense= mat_ref_dense;
sol_upd_coarse = mat_ref_coarse;

index_simp=edge_refined_elem_mapper(M_coarse,M_dense);

% create dense model
img_dense= eidors_obj('image', 'dense image');
img_dense.elem_data= mat_ref_dense;
img_dense.fwd_model= M_dense;
pdense= np_fwd_parameters( M_dense );
vtx_dense      = pdense.vtx;
simp_dense     = pdense.simp;
gnd_ind_dense  = pdense.gnd_ind;
elec_dense     = pdense.elec;
I_dense        = pdense.I;
Ib_dense       = pdense.Ib;
df_dense       = pdense.df;
elec_dense     = pdense.elec;
zc             = pdense.zc;
perm_sym       = pdense.perm_sym;
no_pl          = 1; % FIXME
el_no          = pdense.n_elec;

% create coarse model
img_coarse= eidors_obj('image', 'coarse image');
img_coarse.elem_data= mat_ref_coarse;
img_coarse.fwd_model= M_coarse;
pcoarse= np_fwd_parameters( M_coarse );
vtx_coarse      = pcoarse.vtx;
simp_coarse     = pcoarse.simp;
gnd_ind_coarse  = pcoarse.gnd_ind;
elec_coarse     = pcoarse.elec;
I_coarse        = pcoarse.I;
Ib_coarse       = pcoarse.Ib;
df_coarse       = pcoarse.df;
elec_coarse     = pcoarse.elec;

%Set the tolerance for the forward solver
for iter= 1:inv_model.parameters.max_iterations;                     
    [E_dense,D_dense,Ela_dense,pp_dense] = ...
        fem_master_full(vtx_dense,simp_dense,sol_upd_dense,gnd_ind_dense,elec_dense,zc,perm_sym);
 
 if iter==1
   %sprintf('Current fields for iteration %d',i)
   [V_dense] = forward_solver(E_dense,I_dense,tol_for,pp_dense);
   [viH,viV,indH_dense,indV,df] = get_3d_meas(elec_dense,vtx_dense,V_dense,Ib_dense,no_pl);
   dfv = df(1:2:end);
   vi = (viH); 
   %sprintf('Measurement fields for iteration %d',i)
   [v_f_dense] = m_3d_fields(vtx_dense,el_no,indH_dense,E_dense,tol_for,gnd_ind_dense);
else
   %sprintf('Current fields for iteration %d',i)
   [V_dense] = forward_solver(E_dense,I_dense,tol_for,pp_dense,V_dense);
   [viH,viV,indH_dense,indV,df] = get_3d_meas(elec_dense,vtx_dense,V_dense,Ib_dense,no_pl);
   dfv = df(1:2:end);
   vi = (viH); 
   %sprintf('Measurement fields for iteration %d',i)
   [v_f_dense] = m_3d_fields(vtx_dense,el_no,indH_dense,E_dense,tol_for,gnd_ind_dense,v_f_dense);
end
    
    
%% DENSE MESH TO COARSE MESH HERE for scaler potentials

[V_coarse] = mapvd(V_dense,index_vtx,elec_dense,vtx_coarse,vtx_dense);

[v_f_coarse] = mapvm(v_f_dense,index_vtx,indH_dense,vtx_coarse,vtx_dense);

[viH,viV,indH_coarse,indV,df_coarse] = get_3d_meas(elec_coarse,vtx_coarse,V_coarse,Ib_coarse,no_pl); % for dfv_coarse in Jacobian calculation

dfv_coarse = df_coarse(1:2:end);

[E_coarse,D_coarse,Ela_coarse,pp_coarse] = fem_master_full(vtx_coarse,simp_coarse,mat_ref_coarse,gnd_ind_coarse,elec_coarse,zc,perm_sym);  % for D in Jacobian calculation

[J_coarse] = jacobian_3d_with_fields(V_coarse,Ela_coarse,D_coarse,I_coarse,elec_coarse,vtx_coarse,simp_coarse,gnd_ind_coarse,sol_upd_coarse,zc,v_f_coarse,dfv_coarse,tol_for,perm_sym);

sol_coarse = (J_coarse.'*J_coarse + tfac^2*RtR)\ (J_coarse.' * (voltage - vi));

% COARSE MESH TO DENSE MESH HERE for conductivity

sol_dense=sol_coarse(c2d_elem);

sol_upd_dense = sol_upd_dense + sol_dense;
sol_upd_coarse = sol_upd_coarse + sol_coarse;

sol_coarse = sol_upd_coarse;
sol_dense = sol_upd_dense;

sol_array(:,iter)=sol_coarse;

sprintf('Error norm at iteration %d is %f',iter,norm(voltage - vi))

end


% create a data structure to return
img.name= 'solved by inv_solve_trunc_iterative';
img.elem_data = sol_array;
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

function [V_coarse] = mapvd(V_dense,index_vtx,elec,vtx_coarse,vtx_dense);

% This function maps the scaler potential field
% from the dense mesh to the coarse
% mesh using verticies extraction.

for j=1:size(elec,1);

    for i=1:size(vtx_coarse,1);

        V_coarse(i,j)=V_dense((index_vtx(i,1)),j);
        
    end
    
end

V_coarse_elec=V_dense((size(vtx_dense,1)+1):(size(V_dense,1)),:);

V_coarse=[V_coarse;V_coarse_elec];

% This function maps the scaler potential field
% from the dense mesh to the coarse
% mesh using verticies extraction.
function [v_f_coarse] = mapvm(v_f_dense,index_vtx,indH_dense,vtx_coarse,vtx_dense);


for j=1:size(indH_dense,1);

    for i=1:size(vtx_coarse,1);

        v_f_coarse(i,j)=v_f_dense((index_vtx(i,1)),j);
        
    end
    
end

v_f_coarse_elec=v_f_dense((size(vtx_dense,1)+1):(size(v_f_dense,1)),:);

v_f_coarse=[v_f_coarse;v_f_coarse_elec];
