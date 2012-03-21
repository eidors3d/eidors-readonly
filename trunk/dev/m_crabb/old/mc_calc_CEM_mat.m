function img=mc_calc_CEM_mat(img)
%Takes a forward model and calculates Az, Aw, Ad for CEM, stoing them in
%mdl.solver

%Get the forward model form the image
mdl=img.fwd_model;

%Get the electrode structure, find number of electrodes
elecstruc=mdl.electrode; nelecs=size(elecstruc,2);

%Get the boundary strucutre, find number of boundaries
boundstruc=mdl.bound; nbounds=size(boundstruc,2);

%Get the node structrue, find number of nodes and problem dim
nodestruc=mdl.nodes; nnodes=size(nodestruc,1); nodedim=size(nodestruc,2);

%Loop over boundarys and calculate local Aw/Az matrices
for ii=1:nbounds
    %List by row of coordinates of on the boundaryNodal coordinates on the boundary
    thisb=nodestruc(boundstruc(ii).nodes,:);
    
    %Find the magnitude Jacobian of the mapping in 2D/3D
    %NB:Scalings are consistent with reference element shape
    if(nodedim==2)
        %Length = |(x2-x1)| (x1,x2 vector of coords)
        diff21=thisb(2,:)-thisb(1,:);
        magjacbound=sqrt(diff21(1)^2+diff21(2)^2);
    elseif(nodedim==3)
        %Area = 0.5*|(x3-x1)x(x3-x2)| (x1,x2,x3 vector of coords)
        diff31=thisb(3,:)-thisb(1,:); diff32=thisb(3,:)-thisb(2,:);
        %Cross product and area (double consistent with reference shape)
        diffprod=cross(diff31,diff32);
        magjacbound=sqrt(diffprod(1)^2+diffprod(2)^2+diffprod(3)^2);
    end
    
    %Find quadrature points and weights for integration over boundaries
    %These rules are consistent with the reference element shape choice
    [weight,xcoord,ycoord]=boundgaussquad(boundstruc(ii).type);

    %Initialise Azlocmat/Awlocmat and find local matrices
    %The 0.5 factor comes from choice of the weights above
    Azmat=0; Awmat=0;
    for kk=1:size(weight,2)
        Azmat = Azmat + 0.5*weight(kk)* ...
            (boundshapefunc(boundstruc(ii).type,xcoord(kk),ycoord(kk)))'* ...
            (boundshapefunc(boundstruc(ii).type,xcoord(kk),ycoord(kk)))*magjacbound;
        Awmat = Awmat + 0.5*weight(kk)* ...
            (boundshapefunc(boundstruc(ii).type,xcoord(kk),ycoord(kk)))*magjacbound;
    end              
    %Multiply by the boundary scaling factor and store in mdl.boundary
    boundstruc(ii).Azmat=Azmat; boundstruc(ii).Awmat=Awmat; 
end

%Store the new boundary structure in the forward model
mdl.bound=boundstruc;

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

%Store the matrices in mdl.solver
mdl.solver.Aw=Aw; mdl.solver.Az=Az; mdl.solver.Ad=Ad;

%Put the forward model back into the image
img.fwd_model=mdl;