function[data] = mc_fwd_solve(fwd_model,img)
%Solve for voltages (nodes/electrodes) for a forward model. 

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Modify the forward model to be of my type
[bound,elem,nodes] = fem_modify(fwd_model);
fwd_model.bound=bound; fwd_model.elem=elem; fwd_model.nodes=nodes;
img.fwd_model=fwd_model; %CHANGE THIS

%Calculate the total stiffness matrix
s_mat = calc_system_mat(fwd_model,img); At=s_mat.E;


%Find electrode stucture and no.of electrodes and initialize vector
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);

%Find stim strucutre and no. stimulations
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 

%Find node structure and find no.nodes 
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 

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
    %Assign correct size unknowns and right hand side matrix
    rhsdata=zeros(nnodes,nstims); 
    nodeunknowns=zeros(nnodes,nstims); 
else
    for i=1:nelecs
        %COMPLETE - Assign electrode at bottom of list
        elecnode(i)=nnodes+i;
    end
    %Assign correct size right hand side matrix
    rhsdata=zeros(nnodes+nelecs,nstims); 
    nodeunknowns=zeros(nnodes+nelecs,nstims); 
end

%Assign currents at correct point in rhs matrix using index vector
for ii=1:nstims
    %The vector of current values for stimulation
    curnode=stimstruc(ii).stim_pattern;
    for i=1:nelecs
        rhsdata(elecnode(i),ii)=curnode(i);
    end
end

%Create index vector and eliminate ground node equation from index
groundnode=fwd_model.gnd_node; idx=1:size(At,1); idx(groundnode)=[];

%Solve the simulated linear system with index
nodeunknowns(idx,:)=forward_solver(At(idx,idx),rhsdata(idx,:));


%Find electrode voltages and store in matrix
%Calculate electrode voltages using index vector elecnode(i)
velec=zeros(nelecs,nstims);
for i=1:nelecs
    %This is the indexed entries in nodeunknowns
    velec(i,:)=nodeunknowns(elecnode(i),:);
end

%Get the measured voltages 
if isfield(fwd_model,'absolute') %Volt sum over boundary is zero
    %Initialize electrode voltage measurements vector (and vectroizing index)
    vmeaselec=zeros(nmeass,1); idx=0;
    for ii=1:nstims
        n_meas=nelecs; %Absolute => no. meas = no. elecs??
        %Subtract mean voltage from data, and add to measured
        vmean=mean(velec(:,ii));
        vmeaselec(idx + (1:n_meas) ) = velec(:,ii)-vmean; 
        idx=idx+n_meas; %Increase counter
    end
else %Volt differences between electrodes
    %Initialize electrode voltage measurements vector (and vectroizing index)
    vmeaselec=zeros(nmeass,1); idx=0;
    for ii=1:nstims
        meas_pat=stimstruc(ii).meas_pattern; %Measurement patterns
        n_meas=size(meas_pat,1); %Number of measures
        vmeaselec(idx + (1:n_meas) ) = meas_pat*velec(:,ii); %Diff data
        idx=idx+n_meas; %Increase counter
    end
end


%Return the electrode voltages in data structure
data.meas= vmeaselec;
data.time= -1; % unknown
data.name= 'solved by mc_fwd_solve';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = nodeunknowns(1:nnodes,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = nodeunknowns;             % all, including CEM nodes
end; end

%And now put this as an eidors object
data= eidors_obj('data',data);

end