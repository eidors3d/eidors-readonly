%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('b2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

img=mk_image(fmdl,1);


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
%Gradient (takes nodal volt to elemen Dvolt) of each greens in source in
%each element
DN1 = zeros(n_nodes,2,n_elems); 

for ii=1:n_nodes
    %xii = fmdl.elem_centre(ii,:);
    xii = fmdl.nodes(ii,:);

%    xii_coord = fmdl.nodessqrt(sum(xii.^2));   
    for jj=1:n_nodes
        yjj = fmdl.nodes(jj,:);
     %   yjj_coord = sqrt(sum(yjj.^2));           
        N1(ii,jj) = calc_neumann_func_freespace(xii,yjj);        
    end      
    
    for jj=1:n_elems
        %yjj = fmdl.nodes(jj,:);        
        yjj = fmdl.elem_centre(jj,:);  
%        yjj_coord = sqrt(sum(yjj.^2));           
        DN1(ii,:,jj) = calc_neumann_grad_func_freespace(xii,yjj);
    end    
end

%Pick which element

%for elem_indices=[4,139,576]
for elem_indices=[1,16,63]

    
img0 = img;
img0.elem_data = DN0(:,1,elem_indices);

img1 = img;
img1.elem_data = DN1(:,1,elem_indices);

img2=img;
img2.elem_data = (DN0(:,1,elem_indices)-DN1(:,1,elem_indices));

%%{
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,elem_indices),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,elem_indices),'b*')

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),(N0(:,elem_indices)-N1(:,elem_indices)),'r*')
%}

figure; 
subplot(121); show_fem(img0); colorbar; 
subplot(122); show_fem(img1); colorbar;

figure; 
show_fem(img2); colorbar; 
%}

figure; 
hist(img2.elem_data,100000)
title(sprintf('Error for Greens function source at (x,y) = (%1.2f,%1.2f)',...
    fmdl.elem_centre(elem_indices,1),fmdl.elem_centre(elem_indices,2)))
end


%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),DU1(:,1,1),'r*')