%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

fmdl.gnd_node=1;

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

img=mk_image(fmdl,1);

% sss =  system_mat_1st_order(img)

%Elements and nodes
n_elems = size(fmdl.elems,1);
n_nodes = size(fmdl.nodes,1);


%Space for the analagous analytic freespace/disc N1/2 and DN1/2
N_FEM=zeros(n_nodes,n_nodes); 
%N2=zeros(n_nodes,n_nodes); 
%N3=zeros(n_nodes,n_nodes); 
DN_FEM = zeros(n_nodes,2,n_elems); 
%DN2 = zeros(n_nodes,2,n_elems); 

bound_nodes = unique(fmdl.boundary);

%ii index - nodes at xii coordinate
for ii=1:n_nodes
    %xii = fmdl.elem_centre(ii,:);
    xii = fmdl.nodes(ii,:);

    %jj index - delta function source at yjj
    for jj=1:n_nodes
        yjj = fmdl.nodes(jj,:);
     %   yjj_coord = sqrt(sum(yjj.^2));           
%        NANA(ii,jj) = calc_neumann_func_freespace(xii,yjj);      
        NANA(ii,jj) = calc_neumann_func_disc(xii,yjj);                 
%       N3(ii,jj) = (1+yjj(1))*(1+yjj(2));                         
%        N3(ii,jj) = yjj(1)^2*yjj(2)^2;                         
%        N3(ii,jj) = yjj(2)^2;                                 
    end      
    
    for jj=1:n_elems
        yjj = fmdl.elem_centre(jj,:);  
%        DNANA(ii,:,jj) = calc_neumann_grad_func_freespace(xii,yjj);
        DNANA(ii,:,jj) = calc_neumann_grad_func_disc(xii,yjj);
    end    
end

%N(x,z) - point source at node loc z i.e. each column is vector of
%DN(x,2,z) = gradient w.r.t z of of N(x,z) maps (:,nodes) to (:,2,elems)
NFEM = calc_neumann_func_nodal(img.fwd_model,img);
D_NFEM =calc_grad_potential_nodal(img,NFEM);
D_NANA =calc_grad_potential_nodal(img,NANA);

%DN3 = calc_grad_potential_nodal(img,N3);


figure;
subplot(221); plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NFEM(2,:),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NANA(2,:),'b*')

subplot(222); plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NFEM(50,:),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NANA(50,:),'b*')

subplot(223); plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NFEM(:,2),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NANA(:,2),'b*')

subplot(224); plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NFEM(:,50),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),NANA(:,50),'b*')

for source_indices=[20,50]
    
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),D_NFEM(:,1,source_indices),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),D_NANA(:,1,source_indices),'b*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),DNANA(:,1,source_indices),'g*')

%{
img0x = img;img0y = img;
img0x.elem_data = squeeze(DN3(source_indices,1,:));
img0y.elem_data = squeeze(DN3(source_indices,2,:));

figure; 
subplot(1,2,1); show_fem(img0x,[1,0,0]);% colorbar; 
subplot(1,2,2); show_fem(img0y,[1,0,0]); %colorbar;
%}

end


%{

nnodes = size(img.fwd_model.nodes,1);
s_mat = calc_system_mat(img);
At = s_mat.E(1:nnodes,1:nnodes);


node_index = 75
rhs_fem = At*N0(:,node_index);

N3p = N3;
rhs_ana = At*N3p(:,node_index)

ind_vec = [2:node_index-1,node_index+1:nnodes]

A11 = At(ind_vec,ind_vec)
A12 = At(ind_vec,[1,node_index])
A21 = A12'
A22 = At([1,node_index],[1,node_index])


N0LH = (A11-A12*(A22\A21))*N0(ind_vec,node_index) 
%N0RH = rhs_fem(ind_vec) - A12*(A22\rhs_fem([1,node_index]))

N3LH = (A11-A12*(A22\A21))*N3p(ind_vec,node_index) 
%N3RH = rhs_fem(ind_vec) - A12*(A22\rhs_fem([1,node_index]))

figure; plot(N0LH,'r'); hold on; plot(N3LH,'b')
figure; plot(N0LH-N3LH,'r')

%plot(N3LH-N3RH)
%figure; plot(rhs_ana,'r'); hold on; plot(rhs_fem,'b')

%}
%%{


%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),DU1(:,1,1),'r*')
%}