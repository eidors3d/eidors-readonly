function[data,nodeunknowns] = mc_fwd_solve_point(img)
%Solve for voltages (nodes/electrodes) for a forward model. 

%Get the forward model
mdl=img.fwd_model;

%Find electrode stucture and no.of electrodes and initialize vector
elecstruc=mdl.electrode; nelecs=size(elecstruc,2);

%Find stim strucutre and no. stimulations
stimstruc=mdl.stimulation; nstims=size(stimstruc,2); 

%Find node structure and find no.nodes 
nodestruc=mdl.nodes; nnodes=size(nodestruc,1); 

%Find ground node, system matrix and eliminate ground node equation
groundnode=mdl.gnd_node; At=mdl.solver.Am;
At(groundnode,:)=[]; At(:,groundnode)=[]; 

%Index vector of node numbers of electrodes on boundary
for i=1:nelecs
    elecnode(i)=elecstruc(i).nodes;
end

%Initialize matrix nodeunknowns ((nnodes+nelecs)*nstims)
nodeunknowns=zeros(nnodes,nstims); 

if 1; tic; %NEW QUICKER CODE
%Create empty set of rhs vectors and assign currents on right nodes
rhsdata=zeros(nnodes,nstims);
for ii=1:nstims
    %The vector of current values for stimulation
    curnode=stimstruc(ii).stim_pattern;
    for i=1:nelecs
        rhsdata(elecnode(i),ii)=curnode(i);
    end
end

%Eliminate the ground node equation
rhsdata(groundnode,:)=[];

%Solve all the linear systems using matlab backslash operator
nodetemp=At\rhsdata;

%Store all these in a new matrix
for ii=1:nstims
    nodeunknowns(:,ii)=[nodetemp(1:(groundnode-1),ii); 0; nodetemp((groundnode:end),ii)];
end; fprintf(1,'Forward Solution '); toc;
end     

%Find the voltages on the electrodes (difference data)
totalmeas=0; %Counter
for kk=1:nstims
    %No. measurements for this pattern
    temp=size(stimstruc(kk).meas_pattern,1);
    totalmeas=totalmeas + temp; %Increase counter
end

%Calculate electrode voltages for point (matrix nelecs*nstim)
velec=zeros(nelecs,nstims);
for i=1:nelecs
    %This is the last entries in nodeunknowns
    velec(i,:)=nodeunknowns(elecnode(i),:);
end

%Initialize electrode voltage measurements vector (and vectroizing index)
vmeaselec=zeros(totalmeas,1); idx=0;
for ii=1:nstims %Loop over measurements
    %Measurement matrix
    meas_pat=stimstruc(ii).meas_pattern;
    %No. measurements
    n_meas=size(meas_pat,1);

    %The measured voltages
    vmeaselec(idx + (1:n_meas) ) = meas_pat*velec(:,ii);
    %Increase the storage counter
    idx=idx+n_meas;
end

% create a data structure to return
data.meas= vmeaselec;
data.time= -1; % unknown
data.name= 'solved by mc_fwd_solve';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = v;                % all, including CEM nodes
end; end

%And now put this as an eidors object
data= eidors_obj('data',data);

end