%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('c2C',16);
fmdl = imdl.fwd_model; %Extract model
fmdl = fix_model(fmdl);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

img=mk_image(fmdl,1);

N0 = calc_neumann_function(img.fwd_model,img);
DN0 =calc_grad_potential(img,N0);

n_elems = size(fmdl.elems,1);
n_nodes = size(fmdl.nodes,1);
N1=zeros(n_nodes,n_elems); %greens function with a source in each element
DN1 = zeros(n_elems,2,n_elems); %gradient (takes nodal volt to elemen Dvolt) of each greens in source in each element

for ii=1:n_elems
    xii = fmdl.elem_centre(ii,:);
%    xii_coord = fmdl.nodessqrt(sum(xii.^2));   
    for jjj=1:n_nodes
        yjjj = fmdl.nodes(jjj,:);
     %   yjj_coord = sqrt(sum(yjj.^2));           
        N1(jjj,ii) = calc_N(xii,yjjj);        
    end      
    
    for jj=1:n_elems
        yjj = fmdl.elem_centre(jj,:);  
%        yjj_coord = sqrt(sum(yjj.^2));           
        DN1(jj,:,ii) = calc_grad_N(xii,yjj);
    end    
end

img0 = img;
img0.elem_data = DN0(:,1,4);

img1 = img;
img1.elem_data = DN1(:,1,4);

img2=img
img2.elem_data = (DN0(:,1,4)-DN1(:,1,4))./DN0(:,1,4);

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N0(:,4),'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),N1(:,4),'b*')

figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),(N0(:,4)-N1(:,4)),'r*')


figure; 
subplot(121); show_fem(img0); colorbar; 
subplot(122); show_fem(img1); colorbar;

figure; 
show_fem(img2); colorbar; 
%hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),DU1(:,1,1),'r*')