function [img]=mc_calc_stiffness_mat(img)
%Stiffness matrix, including piecewise conductivity, for EIT

%Get the forward model from the image
mdl=img.fwd_model;

%Cache node structure and find no. of spatial dimensions and nodes
nodestruc=mdl.nodes; nodedim=size(nodestruc,2); nnodes=size(nodestruc,1); 

%Cache element structure and find no. of elements
elemstruc=mdl.elem; nelems=size(elemstruc,2);

for i=1:nelems
    %Find the list of node numbers for each element
    eleminodelist=elemstruc(i).nodes;
    
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
    magjacelem=det(jacobianelem);
        
    %Find quadrature points and weights for integration ove element
    %These rules are consitent with shape choice of reference element
    [weight,xcoord,ycoord,zcoord]=elemgaussquad(elemstruc(i).type);
    
    %Initialise and find elemental stiffness matrices 
    %The 0.5 factor comes from choice of the weights above
    Ammat=0;
    for kk=1:size(weight,2)
        Ammat = Ammat + 0.5*weight(kk)* ...
            (jacobianelem\delemshapefunc(elemstruc(i).type,xcoord(kk),ycoord(kk),zcoord(kk)))'* ...
            (jacobianelem\delemshapefunc(elemstruc(i).type,xcoord(kk),ycoord(kk),zcoord(kk)))*magjacelem;
    end
    %This is element stiffness matrix (and multiply by its conductivity)
    elemstruc(i).stiff=Ammat*img.elem_data(i);
end

%Put the structure back into images forward model
img.fwd_model.elem=elemstruc;

%Initialise stiffness matrix
Agal=zeros(nnodes,nnodes);

%Assemble global stiffness matrix (Silvester's book!!)
for k=1:nelems %loop over elements
    nnodeselems=size(elemstruc(k).nodes,2); %No. of nodes per element
    for i=1:nnodeselems %loop over nodes in element
        for j=1:nnodeselems %loop over nodes in element
            Agal(elemstruc(k).nodes(i),elemstruc(k).nodes(j)) = ...
            Agal(elemstruc(k).nodes(i),elemstruc(k).nodes(j)) + ...
            (elemstruc(k).stiff(i,j));
        end
    end
end

%Store global stiffness matrix in mdl.solver.Am 
img.fwd_model.solver.Am=Agal;
