function [elemstiff] = mc_calc_elem_stiff(fwd_model)
%Elemental stiffness matrices. Each matrix is (essentially) the derivative
%of the forward map with respect to piecewise constant conductivity. This
%is used to compute the Jacobian.

%Cache node structure and find no. of spatial dimensions and nodes
nodestruc=fwd_model.nodes; nodedim=size(nodestruc,2); 

%Cache element structure and find no. of elements
elemstruc=fwd_model.elem; nelems=size(elemstruc,2);

%Find fem type and find quadrature points/weights for integration over
%element consistent with geometry of reference element
eletype=fwd_model.mc_type; 
[weight,xcoord,ycoord,zcoord]=elemgaussquad(eletype);

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
    %This is element stiffness matrix (multiply by conductivity????)
    elemstiff(i).stiff=Ammat; %elemstiff(i).stiff=Ammat*img.elem_data(i); 
end