%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

fmdl.gnd_node=50;

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

img=mk_image(fmdl,1);

% sss =  system_mat_1st_order(img)

%Elements and nodes
n_elems = size(fmdl.elems,1);
n_nodes = size(fmdl.nodes,1);

%N(x,z) - point source at node loc z i.e. each column is vector of
%DN(x,2,z) = gradient w.r.t z of of N(x,z) maps (:,nodes) to (:,2,elems)
N0 = calc_neumann_func_nodal(img.fwd_model,img);
DN0 =calc_grad_potential_nodal(img,N0);

%Space for the analagous analytic freespace/disc N1/2 and DN1/2
N1=zeros(n_nodes,n_nodes); 
N2=zeros(n_nodes,n_nodes); 
%N3=zeros(n_nodes,n_nodes); 
DN1 = zeros(n_nodes,2,n_elems); 
DN2 = zeros(n_nodes,2,n_elems); 

bound_nodes = unique(fmdl.boundary);

%ii index - nodes at xii coordinate
for ii=1:n_nodes
    %xii = fmdl.elem_centre(ii,:);
    xii = fmdl.nodes(ii,:);

    %jj index - delta function source at yjj
    for jj=1:n_nodes
        yjj = fmdl.nodes(jj,:);
     %   yjj_coord = sqrt(sum(yjj.^2));           
        N1(ii,jj) = calc_neumann_func_freespace(xii,yjj);      
        N2(ii,jj) = calc_neumann_func_disc(xii,yjj);          
    end      
    
    for jj=1:n_elems
        %yjj = fmdl.nodes(jj,:);        
        yjj = fmdl.elem_centre(jj,:);  
%        yjj_coord = sqrt(sum(yjj.^2));           
        DN1(ii,:,jj) = calc_neumann_grad_func_freespace(xii,yjj);
        DN2(ii,:,jj) = calc_neumann_grad_func_disc(xii,yjj);
    end    
end

%N0 = calc_neumann_func_nodal(img.fwd_model,img);
%DN2 =calc_grad_potential_nodal(img,N2);


%Put same ground node on analytic solution
%N0=N0'

for iii=1:n_nodes
   N2(:,iii) = N2(:,iii) - N2(fmdl.gnd_node,iii);
end

%sum(N0(bound_nodes,80))

%for node_indices=[5,10,20,50,100]
%for node_indices=[10,15,20,35,50]
for node_indices=[138]

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,node_indices),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,node_indices),'g*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N2(:,node_indices),'b*')

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,node_indices)-N2(:,node_indices),'r*')
%figure;
%plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,node_indices)-N2(:,node_indices),'g*')
%figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N2(:,node_indices)-N3(:,node_indices),'b*')
%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,node_indices)./N2(:,node_indices),'b*')

%figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),...
%    (N0(:,potential_indices)-N1(:,potential_indices)),'r*')    
    
img0x = img;img0y = img;
img0x.elem_data = squeeze(DN0(node_indices,1,:));
img0y.elem_data = squeeze(DN0(node_indices,2,:));

img1x = img;img1y = img;
img1x.elem_data = squeeze(DN1(node_indices,1,:));
img1y.elem_data = squeeze(DN1(node_indices,2,:));

img2x = img;img2y = img;
img2x.elem_data = squeeze(DN2(node_indices,1,:));
img2y.elem_data = squeeze(DN2(node_indices,2,:));

img20x = img;img20y = img;
img20x.elem_data = squeeze(DN2(node_indices,1,:))-squeeze(DN0(node_indices,1,:));
img20y.elem_data = squeeze(DN2(node_indices,2,:))-squeeze(DN0(node_indices,2,:));

img21x = img;img21y = img;
img21x.elem_data = squeeze(DN2(node_indices,1,:))-squeeze(DN1(node_indices,1,:));
img21y.elem_data = squeeze(DN2(node_indices,2,:))-squeeze(DN1(node_indices,2,:));

%img3x=img; img3y=img;
%img3x.elem_data = squeeze((DN0(potential_indices,1,:)-DN1(potential_indices,1,:))...
%    ./abs(DN0(potential_indices,1,:)));
%img3y.elem_data = squeeze((DN0(potential_indices,2,:)-DN2(potential_indices,2,:))...
%    ./abs(DN0(potential_indices,2,:)));

figure; 
subplot(5,2,1); show_fem(img0x,[1,0,0]);% colorbar; 
subplot(5,2,2); show_fem(img0y,[1,0,0]); %colorbar;
subplot(5,2,3); show_fem(img1x,[1,0,0]);% colorbar; 
subplot(5,2,4); show_fem(img1y,[1,0,0]); %colorbar;
subplot(5,2,5); show_fem(img2x,[1,0,0]); %colorbar; 
subplot(5,2,6); show_fem(img2y,[1,0,0]); %colorbar; 
subplot(5,2,7); show_fem(img20x,[1,0,0]); %colorbar; 
subplot(5,2,8); show_fem(img20y,[1,0,0]); %colorbar; 
subplot(5,2,9); show_fem(img21x,[1,0,0]); %colorbar; 
subplot(5,2,10); show_fem(img21y,[1,0,0]); %colorbar; 
%subplot(427); show_fem(img3x,[1,0,0]); %colorbar; 
%subplot(428); show_fem(img3y,[1,0,0]); %colorbar; 

%figure; 
%subplot(121);hist(img2x.elem_data,10000)
%title(sprintf('Error for gradx Greens function source at (x,y) = (%1.2f,%1.2f)',...
%    fmdl.elem_centre(potential_indices,1),fmdl.elem_centre(potential_indices,2)))
%subplot(122);hist(img2y.elem_data,10000)
%title(sprintf('Error for grady Greens function source at (x,y) = (%1.2f,%1.2f)',...
%    fmdl.elem_centre(potential_indices,1),fmdl.elem_centre(potential_indices,2)))
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