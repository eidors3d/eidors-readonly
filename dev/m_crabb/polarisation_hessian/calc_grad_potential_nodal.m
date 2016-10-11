function [ DU0 ] = calc_grad_potential_nodal( img, u0 )
%CALC_GRAD_U0 
% Calculates gradient of npotentials u (defined on nodes) to N gradient Du
% (defined on elements)
% INPUT
% img
% DIM(u0)  nnodes x nodes u0(x,z) = -delta_z dim(u0(x,z)) = (x-nodes,z-nodes)
% OUTPUT
% DIM(DU0) nnodes x dim x nelems Dzu0(x,z) - -Dzu0(x,z) dim(Dzu0(x,z)) = (x-nodes,2,z-elems)

% Calculate gradient of background wavefield

% TODO: needs to be cached really...

fwd_model= img.fwd_model;


%
eletype=fwd_model.approx_type; 
if(strcmp(eletype,'tri3'))
    dim=2; order=0;
elseif(strcmp(eletype,'tri6'))
    dim=2; order=2;
elseif(strcmp(eletype,'tri10'))
    dim=2; order=4;
elseif(strcmp(eletype,'tet4'))
    dim=3; order=0;
elseif(strcmp(eletype,'tet10'))
    dim=3; order=2;
else  
    error('Element type not recognised for integration rules');
end

%Find derivative of shape function on domain element centre, local
%coordintes
if dim==2
    dphi= element_d_shape_function(eletype,1/3.,1/3.,0);
else
    error('dim > 2 not available')
end



%Cache node structure and find no. of spatial dimensions and nodes
%Cache element structure and find no. of elements
nodestruc=fwd_model.nodes; nodedim=size(nodestruc,2); nnodes=size(nodestruc,1); 
elemstruc=fwd_model.elems; nelems=size(elemstruc,1);
n_potentials = size(u0,1);

% Initialise DU0
DU0 = zeros(nnodes, dim, nelems);

%Loop over the elements and calculate local Am matrix
for ii=1:nelems
%     s:-1:1
    %Find the list of node numbers for each element
    eleminodelist=elemstruc(ii,:);
    
    %List by row of coordinate on the element
    thise = nodestruc(eleminodelist,:);
    
    %Find the Jacobian of the mapping in 2D and 3D
    if(nodedim==2); jacobianelem = ... %2D Jacobian of mapping
            [thise(2,1)-thise(1,1),thise(2,2)-thise(1,2); ...
            thise(3,1)-thise(1,1),thise(3,2)-thise(1,2)];  
    elseif(nodedim==3); jacobianelem = ... %3D Jacobian of mapping
            [thise(2,1)-thise(1,1),thise(2,2)-thise(1,2),thise(2,3)-thise(1,3); ...
            thise(3,1)-thise(1,1),thise(3,2)-thise(1,2),thise(3,3)-thise(1,3); ...
            thise(4,1)-thise(1,1),thise(4,2)-thise(1,2),thise(4,3)-thise(1,3)];
    end
    
    %Find the magnitude of the Jacobian of the mapping
    % magjacelem=det(jacobianelem);
    magjacelem=abs(det(jacobianelem));
    
    % global d_shape_function
    dphi_ii = (jacobianelem\dphi).';%*magjacelem;
    
    % Loop over drive patterns
    for jj=1:n_potentials
       DU0(jj,:,ii) = sum(dphi_ii.*repmat(u0(eleminodelist,jj),1,dim), 1 );        
    end        
           
%     %Initialise and find elemental stiffness matrices 
%     Ammat=0;
%     for kk=1:size(weight,2)
%         Ammat = Ammat + weight(kk)* ...
%             (jacobianelem\dphi(:,:,kk))'* ...
%             (jacobianelem\dphi(:,:,kk))*magjacelem;
%     end

%     %SPEED UP
%     %Can we get system_mat_fields here to speed Jacobian?
%     
%     %This is element stiffness matrix (and multiply by its conductivity)
%     stiff=Ammat*img.elem_data(i); 
%     
%     %Assemble global stiffness matrix (Silvester's book!!)    
%     Agal(elemstruc(i,:), elemstruc(i,:)) = Agal(elemstruc(i,:), elemstruc(i,:)) + stiff;

end



end

