function data= np_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% $Id: np_inv_solve.m,v 1.1 2004-07-10 02:40:22 aadler Exp $

% calculate parameters from input structures
keyboard
fwd_model= inv_model.fwd_model;
vtx= fwd_model.nodes;
simp= fwd_model.elems;
elems= size(simp,1);
gnd_ind= fwd_model.gnd_node;

elec= zeros(length(fwd_model.electrode ), ...
            length(fwd_model.electrode(1).nodes) );
zc  = zeros(length(fwd_model.electrode ), 1);

for i=1:length(fwd_model.electrode);
    elec(i,:)= fwd_model.electrode(i).nodes;
    zc(i)    = fwd_model.electrode(i).z_contact;
end


mat_ref= ones(elems,1); % homogeneous background for jacobian

%Set the tolerance for the forward solver
tol = 1e-5;
 [v_f] = m_3d_fields(vtx,32,indH,Eref,tol,gnd_ind);

% Calculating the Jacobian
[J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,dfr,tol,sym);

% Calculating a smoothing prior
[Reg] = iso_f_smooth(simp,vtx,3,1);

% Calculating a linear inverse solution

tfac= inv_model.hyperparameter;
sol = (J'*J +  tfac*Reg'*Reg)\J' * dvaG;


% create a data structure to return
image.name= 'solved by np_inv_solve';
image.elem_data = sol;
image.type = 'real conductivity differences';
image.fwd_model= fwd_model;
image.inv_model= inv_model;
