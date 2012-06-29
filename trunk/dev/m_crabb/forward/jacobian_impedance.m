function J=jacobian_impedance(fwd_model,img)
%Find the Jacobian associated with an image (and forward model) 
%for the contact impedance
%M Crabb - 29.06.2012

%If function called only with image extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Modify the forward model to be of my type
%%fwd_model = mc_fem_modify(fwd_model); img.fwd_model=fwd_model;
[bound,elem,nodes] = fem_modify(fwd_model); 
fwd_model.bound=bound; fwd_model.elem=elem; fwd_model.nodes=nodes;
img.fwd_model=fwd_model;

%Calculate the total stiffness matrix
s_mat = calc_system_mat(fwd_model,img); At=s_mat.E;
 
%Find electrode stucture and no.of electrodes 
%Find stim strucutre and no. stimulations
%Find node structure and find no.nodes 
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 

%Find total number of measurements
nmeass=0;
for k=1:nstims
    stimkmeasmatrix = stimstruc(k).meas_pattern;
    nmeass=nmeass+size(stimkmeasmatrix,1);
end

%Initialise node to electrode matrix
Node2Elec=sparse(nelecs,nnodes+nelecs);
for i=1:nelecs
    %Assign electrode at bottom of list
    elecnode(i)=nnodes+i;
    Node2Elec(i,elecnode(i))=1;
end
        
%Assign correct size unknowns and right hand side matrix (forward)
datafwd=zeros(nnodes+nelecs,nstims); 
nodeunknownsfwd=zeros(nnodes+nelecs,nstims); 

%Loop over stimulations and assign current matrix
for ii=1:nstims
    %The vector of current values for stimulation
    curnode=stimstruc(ii).stim_pattern;
    for i=1:nelecs
        datafwd(elecnode(i),ii)=curnode(i);
    end
end

%Create index vector and eliminate ground node equation from index
groundnode=fwd_model.gnd_node; idx=1:size(At,1); idx(groundnode)=[];

%Solve the simulated linear system with index
nodeunknownsfwd(idx,:)=left_divide(At(idx,idx),datafwd(idx,:));

%Calculate Jacobian tensor - DE_{i,j,k} == dV_i,j / dz_k
%V_i,j - voltage change on electrode i for stim j
%z_k - impedance change on electrode k
DE= zeros(nelecs,nstims,nelecs);

%First step, we only want to pick off the ith electrode
zi2E(:,idx) = Node2Elec(:,idx)/At(idx,idx);

%Calculate the partial derivative matrix for kth change
for k=1:nelecs           
    %Create the FEM derivative matrix
    dA_dzk=mc_calc_d_system_mat(fwd_model,k);

    %Now form product with solution
    DE(:,:,k) = zi2E(:,idx)*dA_dzk(idx,idx)*nodeunknownsfwd(idx,:);
end

%Calculate Jacobian matrix (measurement patterns specified here)
cntjac=0; J=zeros(nmeass,nelecs);
for j=1:nstims   
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   DEj = reshape( DE(:,j,:), nelecs, nelecs);
   J( cntjac+(1:n_meas),: ) = meas_pat*DEj;
   cntjac = cntjac + n_meas;
end; 

%Negative Jacobian for injected currents??
J= -J;  

end

function [At_d]=mc_calc_d_system_mat(fwd_model,elec_no)
%Takes a forward model and calculates the derivatives with respect to the
%contact impedance of the lth electrode

%Get the electrode structure, find number of electrodes
%Get the boundary strucutre, find number of boundaries
%Get the node structrue, find number of nodes and problem dim
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
boundstruc=fwd_model.bound; nbounds=size(boundstruc,2);
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); nodedim=size(nodestruc,2);

%Connect boundary/electrode -pPut boundary into old matrix strucutre
for i=1:nbounds
    boundstrucold(i,:)=boundstruc(i).nodes;
end

%Andy Adler's function, outputs rows of mdl.boundary of electrode and get
%the contact impedance
[bdy_idx,bdy_area]=find_electrode_bdy(boundstrucold(:,1:nodedim),nodestruc,elecstruc(elec_no).nodes);
elecimped=elecstruc(elec_no).z_contact;

%Find fem type and find quadrature points/weights for integration over
%boundaries consistent with geometry of reference boundary
eletype=fwd_model.mc_type; 
[weight,xcoord,ycoord]=boundgaussquad(eletype);

%Loop over boundarys and calculate local Aw/Az matrices
for ii=1:size(bdy_idx,1)
    %List by row of coordinates of on the boundary
    thisb=nodestruc(boundstruc(bdy_idx(ii)).nodes,:);
    
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
end

%1. Initialise global Az_d/Aw_d/Ad_d matrices and assemble a la Silvester
Az_d=zeros(nnodes,nnodes); Aw_d=zeros(nnodes,nelecs); Ad_d=zeros(nelecs,nelecs);

%Find total electrode area (absolute values)
elecarea=0;
for i=1:size(bdy_area,2); elecarea = elecarea + abs(bdy_area(i)); end

%Form the matrix Ad_d
Ad_d(elec_no,elec_no)=-elecarea/elecimped^2;
    
%Loop over each boundary of electrode
for kb=1:size(bdy_idx,1) 
    %Node numbers for this boundary
    boundnodes=boundstruc(bdy_idx(kb)).nodes;
    %Loop over the nodes on the boundary and form Az/Aw matrices 
    for i=1:size(boundnodes,2)
        for j=1:size(boundnodes,2)
            %Form Az matrix in inner loop       
            Az_d(boundnodes(i),boundnodes(j)) = Az_d(boundnodes(i),boundnodes(j)) - Azmat(i,j)/elecimped^2;
        end
        %Form Aw matrix in outer loop        
        Aw_d(boundnodes(i),elec_no) = Aw_d(boundnodes(i),elec_no) + Awmat(i)/elecimped^2;
    end
end


%Form the derivative of the stiffness matrix matrix
At_d=zeros(nnodes+nelecs,nnodes+nelecs);
At_d(1:nnodes,1:nnodes) = Az_d;
At_d(1:nnodes,nnodes+1:nnodes+nelecs) = Aw_d;
At_d(nnodes+1:nnodes+nelecs,1:nnodes)=Aw_d';
At_d(nnodes+1:nnodes+nelecs,nnodes+1:nnodes+nelecs)=Ad_d;

end
