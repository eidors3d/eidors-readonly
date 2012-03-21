function [Agal]=mc_calc_stiffness(fwd_model,img)
%Stiffness matrix, including piecewise conductivity, for EIT. The second
%argument is for Jacobian, it gives discretied gradients in element.

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Cache node structure and find no. of spatial dimensions and nodes
nodestruc=fwd_model.nodes; nodedim=size(nodestruc,2); nnodes=size(nodestruc,1); 

%Cache element structure and find no. of elements
elemstruc=fwd_model.elem; nelems=size(elemstruc,2);

%Find fem type and find quadrature points/weights for integration over
%element consistent with geometry of reference element
eletype=fwd_model.mc_type; 
[weight,xcoord,ycoord,zcoord]=elemgaussquad(eletype);

%Loop over the elements and calculate local Am matrix
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
           
    %Initialise and find elemental stiffness matrices 
    Ammat=0;
    for kk=1:size(weight,2)
        Ammat = Ammat + weight(kk)* ...
            (jacobianelem\delemshapefunc(eletype,xcoord(kk),ycoord(kk),zcoord(kk)))'* ...
            (jacobianelem\delemshapefunc(eletype,xcoord(kk),ycoord(kk),zcoord(kk)))*magjacelem;
    end
    %This is element stiffness matrix (and multiply by its conductivity)
    elemstruc(i).stiff=Ammat*img.elem_data(i); 
end


%Initialise global stiffness matrix
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

%OLD STORAGE - NOT REQUIRED ANYMORE
if 0
    %Put the structure back into images forward model
    img.fwd_model.elem=elemstruc;

    %Store global stiffness matrix in fwd_model.solver.Am 
    img.fwd_model.solver.Am=Agal;
end