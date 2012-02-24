function J = mc_calc_jacobian_CEM(fwd_model,img)
%CHANGE THIS NEED TO PUT FMDL HERE TOO...
%Find the Jacobian associated with an image (and forward model)

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Find electrode stucture and no.of electrodes 
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);

%Find stim strucutre and no. stimulations
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 

%Find node structure and find no.nodes 
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 

%Find element structure and create vector of length no. elements
elemstruc=fwd_model.elem; nelems=size(elemstruc,2);

%Find ground node, system matrix and eliminate ground node equation
groundnode=fwd_model.gnd_node; At=fwd_model.solver.At;
At(groundnode,:)=[]; At(:,groundnode)=[]; 


%Step 1:Resolve the forward problem (CHANGE THIS)

if 1; tic; %NEW QUICKER CODE
nodeunknownsfwd=zeros(nnodes+nelecs,nstims);
%Create empty set of rhs vectors and assign currents on right nodes
datafwd=zeros(nnodes+nelecs,nstims);
for ii=1:nstims
    %The vector of current values for stimulation
    curnode=stimstruc(ii).stim_pattern;
    for i=1:nelecs
        datafwd(i+nnodes,ii)=curnode(i);
    end
end
%Eliminate the ground node equation from data
datafwd(groundnode,:)=[];

%Solve all the linear systems using matlab backslash operator
nodetemp=forward_solver(At,datafwd);
%nodetemp=At\datafwd;

%Store all these in a new matrix
for ii=1:nstims
    nodeunknownsfwd(:,ii)=[nodetemp(1:(groundnode-1),ii); 0; nodetemp((groundnode:end),ii)];
end; fprintf(1,'Forward Solution '); toc;
end   


%Step 2: Solve forward problems on ALL measurement pairs

%Step 2(a):Find total number of measurment pairs
if 1; tic; %NEW QUICKER CODE
%Step 1 : Find total number of measurements
totmeas=0;
for k=1:nstims
    %Measurement matrix for kth stim pattern
    stimkmeasmatrix = stimstruc(k).meas_pattern;
    %Add no. of measurements to counter
    totmeas=totmeas+size(stimkmeasmatrix,1);
end

%Initialise RHS matrix (column vectors of rhs vectors)
data=zeros(nnodes+nelecs,totmeas);

%Step 2(b):Put correct currents in the nodes in data matrix
count=0;
for k=1:nstims
    %Find the measurement matrix for particular stimulation
    stimkmeasmatrix = stimstruc(k).meas_pattern; 

    %Loop over the measurement matrix
    for j=1:size(stimkmeasmatrix,1)
        %Find nodes to stimulate on measurement electrode
        stimkmeasjcur = stimkmeasmatrix(j,:);
        
        %Add on to the counter
        count=count+1;

        %Put currents into RHS vector at the correct indexed node number
        for i=1:nelecs
            data(i+nnodes,count)=stimkmeasjcur(i);
        end
    end
end

%Eliminate the ground node equation from data
data(groundnode,:)=[];

%Solve the sets of linear system (each column is solution of measure drive)
nodetemp=At\data;

%Step 2(c) : Form Jacobian using the elemental inner product u^{e}A^{e}v^{e}
%Initialize counter and Jacobian
cntjac=0; J=zeros(totmeas,nelems);
for k=1:nstims
    %Grab the stimulation nodal unknowns for this stimulation
    stimkunknowns = nodeunknownsfwd(:,k);

    %Loop over the measurement matrix
    for j=1:size(stimkmeasmatrix,1)
        %Add on to the counter
        cntjac=cntjac+1;
 
        %Put ground node (u=0) back into nodeunknowns list
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
            elemsens(ii)=stimknodecpy'*elemstruc(ii).stiff*stimkmeasjnodecpy;
        end
        
        %Store the element sensitivity in a matrix
        J(cntjac,:)=elemsens';  
    end   
end; fprintf(1,'Jacobian Calculation '); toc;
end

%Get the Jacobian and normalize measurements (if field exists)
if isfield(fwd_model,'normalize_measurements')
    if fwd_model.normalize_measurements
 %       J= J ./ (fwd_model.solver.voltmeas(:)*ones(1,nelems));
    end
end

%Negative Jacobian for injected currents??
J= -J;  

end