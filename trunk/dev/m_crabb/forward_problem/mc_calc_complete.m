function [Aw,Az,Ad]=mc_calc_complete(fwd_model,img)
%Takes a forward model and calculates Az, Aw, Ad for complete electrode

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end


%PARAMETERS

%Get the electrode structure, find number of electrodes
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);

%Get the boundary strucutre, find number of boundaries
boundstruc=fwd_model.bound; nbounds=size(boundstruc,2);

%Get the node structrue, find number of nodes and problem dim
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); nodedim=size(nodestruc,2);


%Connect boundary/electrode -Put boundary into old matrix strucutre
for i=1:nbounds
    boundstrucold(i,:)=boundstruc(i).nodes;
end
%Loop over electrodes and store in elecstruc the boundary numbers and areas
for ii=1:nelecs
    %Andy Adler's function, outputs rows of mdl.boundary of electrode
    [bdy_idx,bdy_area]=find_electrode_bdy(boundstrucold(:,1:nodedim),nodestruc,elecstruc(ii).nodes);
    
    %Store boundary numbers, and corresponding areas, in mdl.electrode(ii)
    elecstruc(ii).boundidx=bdy_idx; elecstruc(ii).area=bdy_area;  
end

%Find fem type and find quadrature points/weights for integration over
%boundaries consistent with geometry of reference boundary
eletype=fwd_model.mc_type; 
[weight,xcoord,ycoord]=boundgaussquad(eletype);


%Loop over boundarys and calculate local Aw/Az matrices
%TODO - Loop over elecstruc.boundidx, which has boudnary numbers to save time!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for ii=1:nbounds
    %List by row of coordinates of on the boundaryNodal coordinates on the boundary
    thisb=nodestruc(boundstruc(ii).nodes,:);
    
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
        Azmat = Azmat + weight(kk)* ...
            (boundshapefunc(eletype,xcoord(kk),ycoord(kk)))'* ...
            (boundshapefunc(eletype,xcoord(kk),ycoord(kk)))*magjacbound;
        Awmat = Awmat + weight(kk)* ...
            (boundshapefunc(eletype,xcoord(kk),ycoord(kk)))*magjacbound;
    end              
    %Multiply by the boundary scaling factor and store in fwd_model.boundary
    boundstruc(ii).Azmat=Azmat; boundstruc(ii).Awmat=Awmat; 
end

%1. Initialise global Az/Aw/Ad matrices and assemble a la Silvester
Az=zeros(nnodes,nnodes); Aw=zeros(nnodes,nelecs); Ad=zeros(nelecs,nelecs);
%Loop over the electrodes
for ke=1:nelecs
    %Find boundary numbers (and areas) and contact impedance of electrode
    elecbound=elecstruc(ke).boundidx;
    elecimped=elecstruc(ke).z_contact;
    elecboundarea=elecstruc(ke).area;

    %Find total electrode area (absolute values)
    elecarea=0;
    for i=1:size(elecboundarea,2)
        elecarea = elecarea + abs(elecboundarea(i));
    end
    
    %Form the matrix Ad
    Ad(ke,ke)=elecarea/elecimped;
    
    %Loop over each boundary of electrode
    for kb=1:size(elecbound,1) 
        %Node numbers for this boundary
        boundnodes=boundstruc(elecbound(kb)).nodes;
        %Loop over the nodes on the boundary and form Az/Aw matrices
        for i=1:size(boundnodes,2)
            for j=1:size(boundnodes,2)
                %Form Az matrix in inner loop
                Az(boundnodes(i),boundnodes(j)) = Az(boundnodes(i),boundnodes(j)) + ...
                    boundstruc(elecbound(kb)).Azmat(i,j)/elecimped;
            end
            %Form Aw matrix in outer loop
            Aw(boundnodes(i),ke) = Aw(boundnodes(i),ke) - ...
                boundstruc(elecbound(kb)).Awmat(i);
        end
    end
    %Now divide by the imnpedance of electrode to get Aw
    Aw(:,ke)=Aw(:,ke)/elecimped;
end


%OLD STORAGE - NOT REQUIRED ANYMORE
if 0
    %Store the new boundary structure in the forward model
    fwd_model.bound=boundstruc;
 
    %Store the matrices in fwd_model.solver
    fwd_model.solver.Aw=Aw; fwd_model.solver.Az=Az; fwd_model.solver.Ad=Ad;

    %Put the forward model back into the image
    img.fwd_model=fwd_model;
end
