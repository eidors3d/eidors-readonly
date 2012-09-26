function [s_mat] = calc_system_mat_opt( fwd_model, img )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%[bound,elem,nodes]=fem_1st_to_higher_order(fwd_model);
%bound = fwd_model.bound;
elem = fwd_model.elems;
nodes = fwd_model.nodes;

nodedim=size(nodes,2); nnodes=size(nodes,1); 

nelems=size(elem,2);

Mus = 14;
D = 1.0/(3.0*Mus);
Mua = check_elem_data(fwd_model, img); % 0.15
optical_n = 1.4;
Rx = -1.4399*optical_n^(-2)+0.7099/optical_n+0.6681+0.0636*optical_n;
A = (1+Rx)/(1-Rx);

%Cache node structure and find no. of spatial dimensions and nodes
nodestruc=fwd_model.nodes; nodedim=size(nodestruc,2); nnodes=size(nodestruc,1); 

%Cache element structure and find no. of elements
elemstruc=fwd_model.elems; nelems=size(elemstruc,1);

%Find fem type and find quadrature points/weights for integration over
%element consistent with geometry of reference element
try
    eletype=fwd_model.approx_type;
catch
    eletype='tri3';
end
    
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
for i=1:nelems
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
    magjacelem=abs(det(jacobianelem));
    %Initialise and find elemental stiffness matrices 
    Kmat=0;Cmat=0;
    for kk=1:size(weight1,2)
		Kmat = Kmat + weight1(kk)* ...
            (jacobianelem\dphi(:,:,kk))'* ...
            (jacobianelem\dphi(:,:,kk))*magjacelem;
    end
    for kk=1:size(weight2,2)
        Cmat = Cmat + weight2(kk)* ...
            (phi(:,kk))* ...
            (phi(:,kk))' * magjacelem;
    end
    %This is element stiffness matrix (and multiply by its conductivity)
    stiff=Kmat*D+Cmat*Mua(i); 
    
    %Assemble global stiffness matrix (Silvester's book!!)
    Agal(elemstruc(i,:), elemstruc(i,:)) = Agal(elemstruc(i,:), elemstruc(i,:)) + stiff;
    
    %Store the Cmat without multiplication of Mua for Jacobian
    elemstiff(i).elemstiff=Cmat;
end

if(strcmp(eletype,'tri3'))
    dim=1; order=2;
elseif(strcmp(eletype,'tri6'))
    dim=1; order=4;
elseif(strcmp(eletype,'tri10'))
    dim=1; order=6;
elseif(strcmp(eletype,'tet4'))
    dim=2; order=2;
elseif(strcmp(eletype,'tet10'))
    dim=2; order=4;
else  
    error('Element type not recognised for integration rules');
end
[weight,xcoord,ycoord]=gauss_points(dim,order);
for kk=1:size(weight,2)
    bphi(:,kk) = boundary_shape_function(eletype,xcoord(kk),ycoord(kk))';
end

boundstruc=fwd_model.boundary; nbounds=size(boundstruc,2);

for i=1:nbounds
    %List by row of coordinates of on the boundaryNodal coordinates on the boundary
    thisb=nodestruc(boundstruc(i,:),:);

    %Find the magnitude Jacobian of the mapping in 2D/3D
    %NB:Scalings are consistent with reference element shape
    if(nodedim==2)
        %Jacobian = 0.5*|(x2-x1)| (x1,x2 vector of coords)
        diff21=thisb(2,:)-thisb(1,:);
        magjacbound=0.5*sqrt(diff21(1)^2+diff21(2)^2);
    elseif(nodedim==3)
        %Jacobian = |(x3-x1)x(x3-x2)| (x1,x2,x3 vector of coords)
        diffprod=cross(thisb(3,:)-thisb(1,:),thisb(3,:)-thisb(2,:));
        magjacbound=sqrt(diffprod(1)^2+diffprod(2)^2+diffprod(3)^2);
    end
    
    %Initialise Azlocmat/Awlocmat and find local matrices
    Bmat=0;
    for kk=1:size(weight,2)
        Bmat = Bmat + weight(kk)* ...
            (bphi(:,kk))* ...
            (bphi(:,kk))'*magjacbound;
    end
    
    Agal(boundstruc(i,:),boundstruc(i,:)) = Agal(boundstruc(i,:),boundstruc(i,:)) + Bmat/(2*A);
    
end

s_mat.E=sparse(Agal);

%Store individual stiffness matrices for Jacobian
s_mat.elemstiff=elemstiff;

end

function elem_data = check_elem_data(fwd_model, img);
   elem_data = img.elem_data; 
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('system_mat_1st_order: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(fwd_model, 'coarse2fine');
     c2f = fwd_model.coarse2fine;
     sz_c2f = size(c2f);
     switch sz_elem_data(1)
       case sz_c2f(1); % Ok     
       case sz_c2f(2); elem_data = c2f * elem_data;
       otherwise; error(['system_mat_1st_order: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(fwd_model)
       error(['system_mat_1st_order: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end
end

