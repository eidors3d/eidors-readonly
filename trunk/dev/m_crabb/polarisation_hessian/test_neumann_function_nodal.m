%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

img=mk_image(fmdl,1);

%NE0 = calc_neumann_func(img.fwd_model,img);
%DNE0 =calc_grad_potential(img,N0);

%NO(x,z) for point source at nodal location z i.e. each column is vector of
%nodal potentias for a delta function source at node z
N0 = calc_neumann_func_nodal(img.fwd_model,img);
DN0 =calc_grad_potential_nodal(img,N0);

n_elems = size(fmdl.elems,1);
n_nodes = size(fmdl.nodes,1);
%Greens function with a source in each element
%First argument is potential at nodes and second argument is delta function
%source supported in elements
N1=zeros(n_nodes,n_nodes); 
N2=zeros(n_nodes,n_nodes); 
%Gradient (takes nodal volt to elemen Dvolt) of each greens in source in
%each element
DN1 = zeros(n_nodes,2,n_elems); 
DN2 = zeros(n_nodes,2,n_elems); 

for ii=1:n_nodes
    %xii = fmdl.elem_centre(ii,:);
    xii = fmdl.nodes(ii,:);

%    xii_coord = fmdl.nodessqrt(sum(xii.^2));   
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
        %DN2(ii,:,jj) = calc_neumann_grad_func_freespace(xii,yjj);
    end    
end

for node_indices=[2,5,10,20,50,100]
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,node_indices),'r*')
%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,node_indices),'b*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N2(:,node_indices),'g*')

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,node_indices)-N2(:,node_indices),'r*')
%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,node_indices)./N2(:,node_indices),'b*')

%figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),...
%    (N0(:,potential_indices)-N1(:,potential_indices)),'r*')    
    
end


%Pick which solutions to plot
for potential_indices=[1]

    
img0x = img;img0y = img;
img0x.elem_data = squeeze(DN0(potential_indices,1,:));
img0y.elem_data = squeeze(DN0(potential_indices,2,:));


img1x = img;img1y = img;
img1x.elem_data = squeeze(DN1(potential_indices,1,:));
img1y.elem_data = squeeze(DN1(potential_indices,2,:));
 

img2x=img; img2y=img;
img2x.elem_data = squeeze((DN0(potential_indices,1,:)-DN1(potential_indices,1,:))...
    ./abs(DN0(potential_indices,1,:)));
img2y.elem_data = squeeze((DN0(potential_indices,2,:)-DN1(potential_indices,2,:))...
    ./abs(DN0(potential_indices,2,:)));


%{
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,potential_indices),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,potential_indices),'b*')
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),...
    (N0(:,potential_indices)-N1(:,potential_indices)),'r*')
%}

figure; 
subplot(321); show_fem(img0x,[1,0,0]);% colorbar; 
subplot(322); show_fem(img0y,[1,0,0]); %colorbar;
subplot(323); show_fem(img1x,[1,0,0]);% colorbar; 
subplot(324); show_fem(img1y,[1,0,0]); %colorbar;
subplot(325); show_fem(img2x,[1,0,0]); %colorbar; 
subplot(326); show_fem(img2y,[1,0,0]); %colorbar; 

%}

figure; 
subplot(121);hist(img2x.elem_data,10000)
title(sprintf('Error for gradx Greens function source at (x,y) = (%1.2f,%1.2f)',...
    fmdl.elem_centre(potential_indices,1),fmdl.elem_centre(potential_indices,2)))
subplot(122);hist(img2y.elem_data,10000)
title(sprintf('Error for grady Greens function source at (x,y) = (%1.2f,%1.2f)',...
    fmdl.elem_centre(potential_indices,1),fmdl.elem_centre(potential_indices,2)))
end

%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),DU1(:,1,1),'r*')


