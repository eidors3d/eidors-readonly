function [s_mat]=system_mat_higher_order(fwd_model,img)
%Assemble the total stiffness matrix : s_mat.E=At;
%M Crabb - 29.06.2012
%TODO - Sparse assignment of the matrices
if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling SYSTEM_MAT_HIGHER_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

%Find no. of electrodes and no. of ndoes
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 

%Test - Point/Complete electrodes. Assume no mixed model so test first elect
if(size(elecstruc(1).nodes,2)==1 && size(elecstruc(1).nodes,1)==1) %POINT ELECTRODE
    %IF POINT ELECTRODE
    [At,elemstiff]=mc_calc_stiffness(fwd_model,img);
else %COMPLETE ELECTRODE
    [Am,elemstiff]=mc_calc_stiffness(fwd_model,img);
    
     [Aw,Az,Ad]=mc_calc_complete(fwd_model);
     At=zeros(nnodes+nelecs,nnodes+nelecs);     
%     [i,j,s] = find(Am);
%     At=At+sparse(i,j,s,nnodes+nelecs,nnodes+nelecs);
%     [i,j,s] = find(Az);
%     At=At+sparse(i,j,s,nnodes+nelecs,nnodes+nelecs);
%     [i,j,s] = find(Aw);
%     At=At+sparse(i,j+nnodes,s,nnodes+nelecs,nnodes+nelecs);
%     At=At+sparse(j+nnodes,i,s,nnodes+nelecs,nnodes+nelecs);
%     [i,j,s] = find(Ad);
%     At=At+sparse(i+nnodes,j+nnodes,s,nnodes+nelecs,nnodes+nelecs);
    At(1:nnodes,1:nnodes) = Am+Az;
    At(1:nnodes,nnodes+1:nnodes+nelecs) = Aw;
    At(nnodes+1:nnodes+nelecs,1:nnodes)=Aw';
    At(nnodes+1:nnodes+nelecs,nnodes+1:nnodes+nelecs)=Ad;
end

%Put in structure to be compatibile with eidors
s_mat.E=sparse(At);

%Store individual stiffness matrices for Jacobian
s_mat.elemstiff=elemstiff;
end

%STIFFNESS MATRIX PART
function [Agal,elemstiff]=mc_calc_stiffness(fwd_model,img)
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
[weight,xcoord,ycoord,zcoord]=gauss_points(dim,order);

%Find derivative of shape function on domain element
for kk=size(weight,2):-1:1
    dphi(:,:,kk) = element_d_shape_function(eletype,xcoord(kk),ycoord(kk),zcoord(kk));
end

%Initialise global stiffness matrix
Agal=zeros(nnodes,nnodes); %sparse updating non zero slow

%Initialise structure

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
           
    %Initialise and find elemental stiffness matrices 
    Ammat=0;
    for kk=1:size(weight,2)
%         Ammat = Ammat + weight(kk)* ...
%             (jacobianelem\dphi(:,:,kk))'* ...
%             (jacobianelem\dphi(:,:,kk))*magjacelem;
        jdphitemp=jacobianelem\dphi(:,:,kk);
        Ammat = Ammat + weight(kk)* ...
            (jdphitemp'*jdphitemp)*magjacelem;
    end

    %SPEED UP
    %Can we get system_mat_fields here to speed Jacobian?
    
    %Store the Ammat without multiplication of conductivity for Jacobian
    elemstiff(i).elemstiff=Ammat;
   
    %This is element stiffness matrix (and multiply by its conductivity)

    stiff=Ammat*img.elem_data(i); 
    
    %Assemble global stiffness matrix (Silvester's book!!)    
    Agal(elemstruc(i,:), elemstruc(i,:)) = Agal(elemstruc(i,:), elemstruc(i,:)) + stiff;

end
 
end

%COMPLETE ELECTRODE MATRICES
function [Aw,Az,Ad]=mc_calc_complete(fwd_model)
%Takes a forward model and calculates Az, Aw, Ad for complete electrode

%Get the electrode structure, find number of electrodes
%Get the boundary strucutre, find number of boundaries
%Get the node structrue, find number of nodes and problem dim
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);boundstruc=fwd_model.boundary; nodestruc=fwd_model.nodes; 
nnodes=size(nodestruc,1); nodedim=size(nodestruc,2);

%Connect boundary/electrode -Put boundary into old matrix strucutre
%for i=nbounds:-1:1
%    boundstrucold(i,:)=boundstruc(i).nodes;
%end

%Find quadrature points/weights for integration by switching between cases
eletype=fwd_model.approx_type; 
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

%Find shape function on boundary element
for kk=size(weight,2):-1:1
    phi(:,kk) = boundary_shape_function(eletype,xcoord(kk),ycoord(kk))';
end

%1. Initialise global Az/Aw/Ad matrices and assemble a la Silvester
Az=zeros(nnodes,nnodes); Aw=zeros(nnodes,nelecs); Ad=zeros(nelecs,nelecs); %sparse updating non zero slow

%Loop over the electrodes
for ke=1:nelecs    
    %The boundary numbers and areas, outputs rows of mdl.boundary of electrode
    [bdy_idx,bdy_area]=find_electrode_bdy(boundstruc(:,1:nodedim),nodestruc,elecstruc(ke).nodes);
    
    %Store boundary numbers, and corresponding areas
    boundidx_ke=bdy_idx; area_ke=bdy_area;
    
    %Find contact impedance of electrode
    elecimped=elecstruc(ke).z_contact;   
           
    %Find total electrode area (absolute values)
    elecarea=0;
    for i=1:size(area_ke,2)
        elecarea = elecarea + abs(area_ke(i));
    end
    
    %Form the matrix Ad
    Ad(ke,ke)=elecarea/elecimped; 
    
    
    %Loop over boundarys and calculate Aw/Az matrices
    for ii=1:length(boundidx_ke)
        %List by row of coordinates of on the boundaryNodal coordinates on the boundary
        thisb=nodestruc(boundstruc(boundidx_ke(ii),:),:);
    
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
        Azmat=0; Awmat=0;
        for kk=1:size(weight,2)            
            temphikk = phi(:,kk);
            Azmat = Azmat + weight(kk)* ...
                (temphikk*temphikk')*magjacbound;
            Awmat = Awmat + weight(kk)* ...
                (phi(:,kk))'*magjacbound;            
            %Azmat = Azmat + weight(kk)* ...
            %    (phi(:,kk))* ...
            %    (phi(:,kk))'*magjacbound;
            %Awmat = Awmat + weight(kk)* ...
            %    (phi(:,kk))'*magjacbound;
        end         
        
        %Node numbers for this boundary
        boundnodes=boundstruc(boundidx_ke(ii),:);
        
        %Assemble the matrices
        Az(boundnodes,boundnodes) = Az(boundnodes,boundnodes)+Azmat/elecimped;
        Aw(boundnodes,ke) = Aw(boundnodes,ke) - Awmat'/elecimped;
    end
       
end

end

