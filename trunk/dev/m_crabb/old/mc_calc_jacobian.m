function J = mc_calc_jacobian(fwd_model,img)
%Find the Jacobian associated with an image (and forward model)

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Calculate the total stiffness matrix and elemental stiffness matrices
[At] = mc_calc_system_mat(fwd_model,img);
elemstiff = mc_calc_elem_stiff(fwd_model,img);



%PARAMETERS AND ELECTRODE INDEX VECTORS

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
    datafwd=zeros(nnodes,nstims); nodeunknownsfwd=zeros(nnodes,nstims); 
    
    %Assign correct size right hand size matrix (Jacobian)
    datajac=zeros(nnodes,nmeass);   
else
    for i=1:nelecs
        %COMPLETE - Assign electrode at bottom of list
        elecnode(i)=nnodes+i;
    end
    %Assign correct size unknowns and right hand side matrix (forward)
    datafwd=zeros(nnodes+nelecs,nstims); nodeunknownsfwd=zeros(nnodes+nelecs,nstims); 
    
    %Assign correct size right hand size matrix (Jacobian)
    datajac=zeros(nnodes+nelecs,nmeass);
end

%Find ground node, and eliminate ground node equation
groundnode=fwd_model.gnd_node; At(groundnode,:)=[]; At(:,groundnode)=[]; 



%STEP 1 : RESOLVE FORWARD PROBLEM ON STIM PAIRS



%INITIALISE - Ground nodes and assign currents

%Loop over stimulations and assign currents
for ii=1:nstims
    %The vector of current values for stimulation
    curnode=stimstruc(ii).stim_pattern;
    for i=1:nelecs
        datafwd(elecnode(i),ii)=curnode(i);
    end
end

%Eliminate the ground node equation from data
datafwd(groundnode,:)=[];


%SOLVE - Using backslash/forward_solver and find matrix of unknowns

%tic; 
nodetemp=forward_solver(At,datafwd);%nodetemp=At\rhsdata; 
%fprintf(1,'Forward Solution '); toc;

%Store all these in a new matrix
for ii=1:nstims
    nodeunknownsfwd(:,ii)=[nodetemp(1:(groundnode-1),ii); 0; nodetemp((groundnode:end),ii)];
end   




%STEP 2 : SOLVE FORWARD PROBLEM ON ALL MEASUREMENT PAIRS

%INITIALISE - Ground nodes and assign currents
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

%Eliminate the ground node equation from data
datajac(groundnode,:)=[];



%SOLVE - Using backslash/forward_solver

%tic;
nodetemp=At\datajac; %nodetemp=forward_solver(At,datajac);
%fprintf(1,'Jacobian Forward Solution '); toc;



%POST-PROCESS - Form Jacobian using elemental inner product u^{e}A^{e}v^{e}
%Initialize counter and Jacobian
%tic; 
cntjac=0; J=zeros(nmeass,nelems);
for k=1:nstims
    %Grab the stimulation nodal unknowns for this stimulation
    stimkunknowns = nodeunknownsfwd(:,k);

    %Loop over the measurement matrix
    for j=1:size(stimkmeasmatrix,1)
        
        %Put ground node (u=0) back into nodeunknowns list and ++counter
        cntjac=cntjac+1;
        stimkmeasjunknowns=[nodetemp(1:(groundnode-1),cntjac); 0; nodetemp(groundnode:end,cntjac)];

        %Initialise a sensitvity vector
        elemsens=zeros(nelems,1);

        %Loop over and calculate inner product ve'*Ae*ue
        for ii=1:nelems
            %List the nodes of the element
            eleminodes = elemstruc(ii).nodes;
            %Use nodal index to get unknowns at nodal position in element
            stimknodecpy = stimkunknowns(eleminodes,:);
            stimkmeasjnodecpy = stimkmeasjunknowns(eleminodes,:);
            %Calculate inner product
            elemsens(ii)=stimknodecpy'*elemstiff(ii).stiff*stimkmeasjnodecpy;
        end
        
        %Store the element sensitivity in a matrix
        J(cntjac,:)=elemsens';  
    end   
end; 
%fprintf(1,'Jacobian Calculation '); toc;

%Get the Jacobian and normalize measurements (if field exists)
if isfield(fwd_model,'normalize_measurements')
    data=fwd_solve( img );   
    J= J ./ (data.meas(:)*ones(1,nelems));
end

%Negative Jacobian for injected currents??
J= -J;  

end