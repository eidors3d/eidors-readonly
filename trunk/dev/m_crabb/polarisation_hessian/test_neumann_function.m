%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('c2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

img=mk_image(fmdl,1);

N0 = calc_neumann_func(img.fwd_model,img);
DN0 =calc_grad_potential(img,N0);

n_elems = size(fmdl.elems,1);
n_nodes = size(fmdl.nodes,1);
%Greens function with a source in each element
%First argument is potential at nodes and second argument is delta function
%source supported in elements
N1=zeros(n_nodes,n_elems); 
%Gradient (takes nodal volt to elemen Dvolt) of each greens in source in
%each element
DN1 = zeros(n_elems,2,n_elems); 

for ii=1:n_elems
    xii = fmdl.elem_centre(ii,:);
%    xii_coord = fmdl.nodessqrt(sum(xii.^2));   
    for jjj=1:n_nodes
        yjjj = fmdl.nodes(jjj,:);
     %   yjj_coord = sqrt(sum(yjj.^2));           
        N1(jjj,ii) = calc_neumann_func_freespace(xii,yjjj);        
    end      
    
    for jj=1:n_elems
        yjj = fmdl.elem_centre(jj,:);  
%        yjj_coord = sqrt(sum(yjj.^2));           
        DN1(jj,:,ii) = calc_neumann_grad_func_freespace(xii,yjj);
    end    
end

%Pick which element

for elem_indices=[4,139,576]

img0 = img;
img0.elem_data = DN0(:,1,elem_indices);

img1 = img;
img1.elem_data = DN1(:,1,elem_indices);

img2=img;
img2.elem_data = (DN0(:,1,elem_indices)-DN1(:,1,elem_indices))./DN1(:,1,elem_indices);

%{
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,4),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,4),'b*')

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),(N0(:,4)-N1(:,4)),'r*')


figure; 
subplot(121); show_fem(img0); colorbar; 
subplot(122); show_fem(img1); colorbar;

figure; 
show_fem(img2); colorbar; 
%}

figure; 
hist(img2.elem_data,1000)
title(sprintf('Error for Greens function source at (x,y) = (%1.2f,%1.2f)',...
    fmdl.elem_centre(elem_indices,1),fmdl.elem_centre(elem_indices,2)))
end


%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),DU1(:,1,1),'r*')