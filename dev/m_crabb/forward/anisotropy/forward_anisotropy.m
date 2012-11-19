%Choose a 2D or 3D model
m_dim=2; 

%Make an inverse model and extract forward model
if(m_dim==2)
    imdl = mk_common_model('a2C2',16);
else
    imdl = mk_common_model('n3r2',[16,2]);    
end
fmdl = imdl.fwd_model;

%Assign a matrix to each elem data
n_elem=size(fmdl.elems,1); 
if(m_dim==2)
    for i=1:n_elem 
        %Symmetric conductivity tensor Cartesian coordinates
        elem_data(i,1,1,1) =1;
        elem_data(i,1,1,2) =0; elem_data(i,1,2,1) =0;
        elem_data(i,1,2,2) =2;
    end
else
    for i=1:n_elem 
        %Symmetric conductivity tensor Cartesian coordinates
        elem_data(i,1,1,1) =1;
        elem_data(i,1,1,2) =0; elem_data(i,1,2,1) =0;
        elem_data(i,1,1,3) =0; elem_data(i,1,3,1) =0;
        elem_data(i,1,2,2) =3;
        elem_data(i,1,2,3) =0; elem_data(i,1,3,2) =0;
        elem_data(i,1,3,3) =3;   
    end
end

%High-order EIDORS solver
%Change default eidors solvers
fmdl.system_mat = @system_mat_higher_order_anisotropy;

%Add element type
if(m_dim==2)
    fmdl.approx_type    = 'tri3'; % linear
else
    fmdl.approx_type    = 'tet4'; % linear
end
%Make an image and get voltages using high order solver
img1 = mk_image_anisotropy(fmdl,elem_data);
img1.fwd_solve.get_all_meas = 1; %Internal voltage
v1 = fwd_solve(img1); 
v1e=v1.meas; v1all=v1.volt;

%Plot electrode voltages and difference
figure; plot([v1e]);
legend('1'); 

%Plot the voltage distribution
v1all = v1all; 
img1n = rmfield(img1,'elem_data');
img1n.node_data = v1all(1:size(fmdl.nodes,1),1); %add first stim data
figure; show_fem(img1n,1);

%Plot the current distribution
img_v = img1;
img_v.fwd_model.mdl_slice_mapper.npx = 64;
img_v.fwd_model.mdl_slice_mapper.npy = 64;
img_v.fwd_model.mdl_slice_mapper.level = [inf,inf,1.0];
figure; show_current(img_v, v1.volt(:,1));

%Calculate the Jacobian and compare with simple perturbation
img1.fwd_model.jacobian=@jacobian_adjoint_higher_order_anisotropy;
J = calc_jacobian(img1);
Jp=jacobian_adjoint_higher_order_anisotropy_perturb(img1.fwd_model,img1);
norm(Jp-J)/norm(J)

