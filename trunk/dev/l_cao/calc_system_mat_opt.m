function [s_mat] = calc_system_mat_opt( fwd_model, img )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%[bound,elem,nodes]=fem_1st_to_higher_order(fwd_model);
bound = fwd_model.bound;
elem = fwd_model.elem;
nodes = fwd_model.nodes;

nodedim=size(nodes,2); nnodes=size(nodes,1); 

nelems=size(elem,2);

Mus = 10;
D = 1.0/(3.0*Mus);
Mua = 0.1;
optical_n = 1.4;
Rx = -1.4399*optical_n^(-2)+0.7099/optical_n+0.6681+0.0636*optical_n;
A = (1+Rx)/(1-Rx);

%Cache node structure and find no. of spatial dimensions and nodes
nodestruc=fwd_model.nodes; nodedim=size(nodestruc,2); nnodes=size(nodestruc,1); 

%Cache element structure and find no. of elements
elemstruc=fwd_model.elem; nelems=size(elemstruc,2);

%Find fem type and find quadrature points/weights for integration over
%element consistent with geometry of reference element
eletype=fwd_model.approx_type; 
[weight,xcoord,ycoord,zcoord]=element_gauss_points(eletype);
for kk=1:size(weight,2)
    dphi(:,:,kk) = element_d_shape_function(eletype,xcoord(kk),ycoord(kk),zcoord(kk));
    phi(:,kk) = element_shape_function(eletype,xcoord(kk),ycoord(kk),zcoord(kk));
end

%Initialise global stiffness matrix
Agal=zeros(nnodes,nnodes); %sparse updating non zero slow

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
    Kmat=0;Cmat=0;
    for kk=1:size(weight,2)
        Kmat = Kmat + weight(kk)* ...
            (jacobianelem\dphi(:,:,kk))'* ...
            (jacobianelem\dphi(:,:,kk))*magjacelem;
        Cmat = Cmat + weight(kk)* ...
            (phi(:,kk))'* ...
            (phi(:,kk)) * magjacelem;
    end
    %This is element stiffness matrix (and multiply by its conductivity)
    stiff=Kmat*D+Cmat*Mua;%img.elem_data(i); 
    
    %Assemble global stiffness matrix (Silvester's book!!)
    Agal(elemstruc(i).nodes, elemstruc(i).nodes) = Agal(elemstruc(i).nodes, elemstruc(i).nodes) + stiff;
    
end




s_mat.E=sparse(Agal);


end


% %COMPLETE ELECTRODE MATRICES
% function [Aw,Az,Ad]=mc_calc_complete(fwd_model,img)
% %Takes a forward model and calculates Az, Aw, Ad for complete electrode
% 
% %If function called only with image, extract forward model
% if(nargin==1)
%     img=fwd_model; fwd_model=img.fwd_model;
% end
% 
% %Get the electrode structure, find number of electrodes
% %Get the boundary strucutre, find number of boundaries
% %Get the node structrue, find number of nodes and problem dim
% elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
% boundstruc=fwd_model.bound; nbounds=size(boundstruc,2);
% nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); nodedim=size(nodestruc,2);
% 
% %Connect boundary/electrode -Put boundary into old matrix strucutre
% for i=1:nbounds
%     boundstrucold(i,:)=boundstruc(i).nodes;
% end
% 
% %Find fem type and find quadrature points/weights for integration over
% %boundaries consistent with geometry of reference boundary
% eletype=fwd_model.approx_type; 
% [weight,xcoord,ycoord]=boundary_gauss_points(eletype);
% 
% 
% %1. Initialise global Az/Aw/Ad matrices and assemble a la Silvester
% Az=zeros(nnodes,nnodes); Aw=zeros(nnodes,nelecs); Ad=zeros(nelecs,nelecs); %sparse updating non zero slow
% 
% %Loop over the electrodes
% for ke=1:nelecs    
%     %The boundary numbers and areas, outputs rows of mdl.boundary of electrode
%     [bdy_idx,bdy_area]=find_electrode_bdy(boundstrucold(:,1:nodedim),nodestruc,elecstruc(ke).nodes);
%     
%     %Store boundary numbers, and corresponding areas
%     boundidx_ke=bdy_idx; area_ke=bdy_area;
%     
%     %Find contact impedance of electrode
%     elecimped=elecstruc(ke).z_contact;   
%            
%     %Find total electrode area (absolute values)
%     elecarea=0;
%     for i=1:size(area_ke,2)
%         elecarea = elecarea + abs(area_ke(i));
%     end
%     
%     %Form the matrix Ad
%     Ad(ke,ke)=elecarea/elecimped; 
%     
%     
%     %Loop over boundarys and calculate Aw/Az matrices
%     for ii=1:length(boundidx_ke)
%         %List by row of coordinates of on the boundaryNodal coordinates on the boundary
%         thisb=nodestruc(boundstruc(boundidx_ke(ii)).nodes,:);
%     
%         %Find the magnitude Jacobian of the mapping in 2D/3D
%         %NB:Scalings are consistent with reference element shape
%         if(nodedim==2)
%             %Jacobian = 0.5*|(x2-x1)| (x1,x2 vector of coords)
%             diff21=thisb(2,:)-thisb(1,:);
%             magjacbound=0.5*sqrt(diff21(1)^2+diff21(2)^2);
%         elseif(nodedim==3)
%             %Jacobian = |(x3-x1)x(x3-x2)| (x1,x2,x3 vector of coords)
%             diffprod=cross(thisb(3,:)-thisb(1,:),thisb(3,:)-thisb(2,:));
%             magjacbound=sqrt(diffprod(1)^2+diffprod(2)^2+diffprod(3)^2);
%         end
% 
%         %Initialise Azlocmat/Awlocmat and find local matrices
%         Azmat=0; Awmat=0;
%         for kk=1:size(weight,2)
%             Azmat = Azmat + weight(kk)* ...
%                 (boundary_shape_function(eletype,xcoord(kk),ycoord(kk)))'* ...
%                 (boundary_shape_function(eletype,xcoord(kk),ycoord(kk)))*magjacbound;
%             Awmat = Awmat + weight(kk)* ...
%                 (boundary_shape_function(eletype,xcoord(kk),ycoord(kk)))*magjacbound;
%         end         
%         
%         %Node numbers for this boundary
%         boundnodes=boundstruc(boundidx_ke(ii)).nodes;
%         
%         Az(boundnodes,boundnodes) = Az(boundnodes,boundnodes)+Azmat/elecimped;
%         Aw(boundnodes,ke) = Aw(boundnodes,ke) - Awmat'/elecimped;
%         
% %         dofAzi = zeros(1,size(boundnodes,2)*size(boundnodes,2));
% %         dofAzj = zeros(1,size(boundnodes,2)*size(boundnodes,2));
% %         dofAzv = zeros(1,size(boundnodes,2)*size(boundnodes,2));
% %         
% %         dofAwi = zeros(1,size(boundnodes,2));
% %         dofAwj = zeros(1,size(boundnodes,2));
% %         dofAwv = zeros(1,size(boundnodes,2));
% %         
% %         sparseindex = 1;
% %         for i=1:size(boundnodes,2)
% %             for j=1:size(boundnodes,2)
% %                 dofAzi(sparseindex) = boundnodes(i);
% %                 dofAzj(sparseindex) = boundnodes(j);
% %                 dofAzv(sparseindex) = Azmat(i,j)/elecimped;
% %                 sparseindex = sparseindex + 1;
% %             end
% %             dofAwi(i) = boundnodes(i);
% %             dofAwj(i) = ke;
% %             dofAwv(i) = Awmat(i)/elecimped;
% %         end          
% %         Az = Az + sparse(dofAzi,dofAzj,dofAzv,nnodes,nnodes);      
% %         Aw = Aw + sparse(dofAwi,dofAwj,dofAwz,nnodes,nelecs);
%         
%         %Loop over the nodes on the boundary and form Az/Aw matrices
% %         for i=1:size(boundnodes,2)
% %             for j=1:size(boundnodes,2)
% %                 %Form Az matrix in inner loop
% %                 %Az(boundnodes(i),boundnodes(j)) =
% %                 Az(boundnodes(i),boundnodes(j)) + Azmat(i,j)/elecimped;
% %             end
% %             %Form Aw matrix in outer loop
% %             Aw(boundnodes(i),ke) = Aw(boundnodes(i),ke) - Awmat(i)/elecimped;
% %         end        
% 
%     end
%        
% end

