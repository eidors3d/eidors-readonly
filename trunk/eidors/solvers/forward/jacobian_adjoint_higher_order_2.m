function J = jacobian_adjoint_higher_order_2(fwd_model,img)
%Find the Jacobian associated with an image (and forward model)
%Discretization of derivative method
%M Crabb - 29.06.2012

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Modify the forward model to be of my type
%%fwd_model = mc_fem_modify(fwd_model); img.fwd_model=fwd_model;
[bound,elem,nodes] = fem_modify(fwd_model); 
fwd_model.bound=bound; fwd_model.elem=elem; fwd_model.nodes=nodes;
img.fwd_model=fwd_model; %CHANGE THIS!!!!!!!!!!!!!!

%Calculate the total stiffness matrix and elemental stiffness matrices
s_mat = calc_system_mat(fwd_model,img); At=s_mat.E;
%Get the stiffness matrix structure
elemstiff = mc_calc_elem_stiff(fwd_model); %CHANGE THIS!!!!!!!!!!!!!!


%Find electrode stucture and no.of electrodes 
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);

%Find stim strucutre and no. stimulations
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 

%Find node structure and find no.nodes 
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 

%Find element structure and create vector of length no. elements
elemstruc=fwd_model.elem; nelems=size(elemstruc,2);

%Find total number of measurements
nmeass=0;
for k=1:nstims
    stimkmeasmatrix = stimstruc(k).meas_pattern;
    nmeass=nmeass+size(stimkmeasmatrix,1);
end

%Complete or Point? Check first electrode (no mized models) and this changes
%the index vector of what 'node' number corresponds to an electrode
elecnode=zeros(1,nelecs);
if(size(elecstruc(1).nodes,2)==1 && size(elecstruc(1).nodes,1)==1) %POINT
    for i=1:nelecs
        %POINT - Assign electrode index at correct node
        elecnode(i)=elecstruc(i).nodes;
    end
    %Assign correct size unknowns and right hand side matrix (forward)
    datafwd=zeros(nnodes,nstims); 
    nodeunknownsfwd=zeros(nnodes,nstims); 
    
    %Assign correct size unknowns and right hand size matrix (Jacobian)
    datajac=zeros(nnodes,nmeass);  
    nodeunknownsjac=zeros(nnodes,nmeass);
else
    for i=1:nelecs
        %COMPLETE - Assign electrode at bottom of list
        elecnode(i)=nnodes+i;
    end
    %Assign correct size unknowns and right hand side matrix (forward)
    datafwd=zeros(nnodes+nelecs,nstims); 
    nodeunknownsfwd=zeros(nnodes+nelecs,nstims); 
    
    %Assign correct size unknowns and right hand size matrix (Jacobian)
    datajac=zeros(nnodes+nelecs,nmeass);
    nodeunknownsjac=zeros(nnodes+nelecs,nmeass);
end

%STEP 1 : SOLVE FORWARD PROBLEM ON STIM PAIRS

%Loop over stimulations and assign currents
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


%STEP 2 : SOLVE FORWARD PROBLEM ON ALL MEASUREMENT PAIRS

%Loop over stimulations and assign currents
count=0;
for k=1:nstims
    %Find the measurement matrix for particular stimulation
    stimkmeasmatrix = stimstruc(k).meas_pattern; 

    %For each stimulation loop over measurment matrix
    for j=1:size(stimkmeasmatrix,1)
        %Find nodes to stimulate on measurement electrode and ++counter
        stimkmeasjcur = stimkmeasmatrix(j,:); count=count+1;

        %Put currents into RHS vector at the correct indexed node number
        for i=1:nelecs
            datajac(elecnode(i),count)=stimkmeasjcur(i);
        end
    end
end

%Solve the simulated linear system with index
nodeunknownsjac(idx,:)=left_divide(At(idx,idx),datajac(idx,:));


%STEP 3 : FORM JACOBIAN USING INNER PRODUCT u^{e}A^{e}v^{e}

%Initialise measurement counter and Jacobian
cntjac=0; J=zeros(nmeass,nelems);
for k=1:nstims
    %Find the measurement matrix for particular stimulation
    stimkmeasmatrix = stimstruc(k).meas_pattern; 

    %Loop over the measurement matrix
    for j=1:size(stimkmeasmatrix,1)
        %Increase measurement counter
        cntjac=cntjac+1;

        %Loop over and calculate inner product ve'*Ae*ue
        for ii=1:nelems        
            %Global node numbers of element
            eleminodes = elemstruc(ii).nodes;
            %Calculate inner product for Jacobian
            J(cntjac,ii) = nodeunknownsfwd(eleminodes,k)'*elemstiff(ii).stiff*nodeunknownsjac(eleminodes,cntjac);
        end
    end   
end; 

%Get the Jacobian and normalize measurements (if field exists)
if isfield(fwd_model,'normalize_measurements')
%    data=mc_fwd_solve( img );   
%    J= J ./ (data.meas(:)*ones(1,nelems));
end

%Negative Jacobian for injected currents
J= -J;  

end

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
eletype=fwd_model.approx_type; 
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

end
