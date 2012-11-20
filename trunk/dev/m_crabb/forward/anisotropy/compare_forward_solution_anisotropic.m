%FORWARD SOLUTIONS WITH ISOTROPY AND ANISOTROPY IN 2D
%Choose 
for int_ani=[1 2]; %1 is exterior isotropic, interior aniostropic and 2 vice versa

%Make an inverse model and extract forward model
imdl = mk_common_model('c2C2',16);
fmdl = imdl.fwd_model;

%Assign a matrix to each elem data
n_elem=size(fmdl.elems,1); 
if(int_ani==1)
    for i=1:64 
        %Symmetric conductivity tensor Cartesian coordinates
        elem_data(i,1,1,1) =10;
        elem_data(i,1,1,2) =9; elem_data(i,1,2,1) =9;
        elem_data(i,1,2,2) =10;
    end
    %Change for interior pixels    
    for i=64+1:n_elem %Change exterior to isotropic for 'c2C2' model
        %Symmetric conductivity tensor Cartesian coordinates
        elem_data(i,1,1,1) =1;
        elem_data(i,1,1,2) =0; elem_data(i,1,2,1) =1;
        elem_data(i,1,2,2) =1;
    end        
else
    for i=1:n_elem
        %Symmetric conductivity tensor Cartesian coordinates
        elem_data(i,1,1,1) =10;
        elem_data(i,1,1,2) =9; elem_data(i,1,2,1) =9;
        elem_data(i,1,2,2) =10;
    end
    %Change for interior pixels    
    for i=1:64 %Change interior to isotropic
        %Symmetric conductivity tensor Cartesian coordinates
        elem_data(i,1,1,1) =1;
        elem_data(i,1,1,2) =0; elem_data(i,1,2,1) =1;
        elem_data(i,1,2,2) =1;
    end        
end


%High-order EIDORS solver
%Change default eidors solvers
fmdl.system_mat = @system_mat_higher_order_anisotropy;
fmdl.approx_type    = 'tri3'; % linear

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
end
