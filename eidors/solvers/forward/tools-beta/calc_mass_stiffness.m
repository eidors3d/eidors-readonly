%STIFFNESS MATRIX PART
function [Agal,elemstiff]=mc_calc_stiffness2(fwd_model,img)
%Stiffness matrix, including piecewise conductivity, for EIT. The second
%argument is for Jacobian, it gives discretied gradients in element.

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Cache node structure and find no. of spatial dimensions and nodes
%Cache element structure and find no. of elements
nodestruc=fwd_model.nodes; nodedim=size(nodestruc,2); nnodes=size(nodestruc,1); 
elemstruc=fwd_model.elems; nelems=size(elemstruc,1);

%Find quadrature points/weights for integration by switching between cases
eletype=fwd_model.approx_type; 
if(strcmp(eletype,'tri3'))
    dim=2; order1=0; order2=2;
elseif(strcmp(eletype,'tri6'))    
    dim=2; order1=2; order2=4;
elseif(strcmp(eletype,'tri10'))
    dim=2; order1=4; order2=7;
elseif(strcmp(eletype,'tet4'))
    dim=3; order1=0; order2=2;
elseif(strcmp(eletype,'tet10'))
    dim=3; order1=2; order2=4;
else  
    error('Element type not recognised for integration rules');
end
[weight1,xcoord1,ycoord1,zcoord1]=gauss_points(dim,order1);
for kk=1:size(weight1,2)
    dphi(:,:,kk) = element_d_shape_function(eletype,xcoord1(kk),ycoord1(kk),zcoord1(kk));
end

[weight2,xcoord2,ycoord2,zcoord2]=gauss_points(dim,order2);
for kk=1:size(weight2,2)
    phi(:,kk) = element_shape_function(eletype,xcoord2(kk),ycoord2(kk),zcoord2(kk));
end

%Initialise global stiffness matrix
Agal=zeros(nnodes,nnodes); %sparse updating non zero slow

%Loop over the elements and calculate local Am matrix
for i=nelems:-1:1
    %Find the list of node numbers for each element
    eleminodelist=elemstruc(i,:);
    
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
           
    %Elemental stiffness matrices 
    Ammat=0;
    for kk=1:size(weight1,2)
        Ammat = Ammat + weight1(kk)* ...
            (jacobianelem\dphi(:,:,kk))'* ...
            (jacobianelem\dphi(:,:,kk))*magjacelem;
    end

    %Element mass matrices
    Mmmat=0;
    for kk=1:size(weight2,2)
        Mmmat = Mmmat + weight2(kk)* ...
            (phi(:,kk))* ...
            (phi(:,kk))' * magjacelem;
    end
    
    %Store the Ammat without multiplication of conductivity for Jacobian
    elemstiff(i).elemstiff = Ammat;
    elemstiff(i).elemmass  = Mmmat;
   
    %This is element stiffness matrix (and multiply by its conductivity)
    stiff=Ammat*img.elem_data(i); 
    
    %Assemble global stiffness matrix (Silvester's book!!)    
    Agal(elemstruc(i,:), elemstruc(i,:)) = Agal(elemstruc(i,:), elemstruc(i,:)) + stiff;

end

end