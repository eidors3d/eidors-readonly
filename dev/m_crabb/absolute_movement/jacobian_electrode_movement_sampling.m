function [J] = jacobian_electrode_movement_sampling(fwd_model,img)
%This only will work in 2D currently!!!

%Test the arguments
if nargin == 1
   img= fwd_model;
else
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_ADJOINT_HIGHER_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

%Modify the forward model to be of my type
if ~isfield(fwd_model,'approx_type')    || ...
   strcmp(fwd_model.approx_type,'tri3') || ...
   strcmp(fwd_model.approx_type,'tet4')
    %Do nothing
else
    [bound,elem,nodes] = fem_1st_to_higher_order(fwd_model);
    fwd_model.boundary=bound; fwd_model.elems=elem; fwd_model.nodes=nodes;
    img.fwd_model=fwd_model; %CHANGE THIS
end

%Calculate the total stiffness matrix and elemental stiffness matrices
s_mat = calc_system_mat(img); At=s_mat.E;

%Get all the nodes, elems, elecs and boundary info
elecstruc = img.fwd_model.electrode; n_elec = length(elecstruc);
elemstruc = img.fwd_model.elems;     n_elem=length(elemstruc(:,1));
nodestruc =  img.fwd_model.nodes;    nodedim = length(nodestruc(1,:)); n_node = length(nodestruc(:,1));
boundstruc = img.fwd_model.boundary; n_bound = length(boundstruc(:,1));
stimstruc = img.fwd_model.stimulation; n_stim=length(stimstruc);

%Find ground node, system matrix and eliminate ground node equation
groundnode=fwd_model.gnd_node; At(groundnode,:)=[]; At(:,groundnode)=[]; 

%Resolve the forward problem for each stimulation
nodeunknownsfwd=zeros(n_node+n_elec,n_stim);
datafwd=zeros(n_node+n_elec,n_stim);
totmeas=0;
for ii=1:n_stim
    %The vector of current values for stimulation
    curnode=stimstruc(ii).stim_pattern;
    for i=1:n_elec
        datafwd(i+n_node,ii)=curnode(i);
    end
    %Count total measurement number
    stimkmeasmatrix = stimstruc(ii).meas_pattern;
    %Add no. of measurements to counter
    totmeas=totmeas+size(stimkmeasmatrix,1);
end
%Eliminate the ground node equation from data
datafwd(groundnode,:)=[];
%Solve all the linear systems using matlab backslash operator
nodetemp=left_divide(At,datafwd);
%Store all these in a new matrix
for ii=1:n_stim
    nodeunknownsfwd(:,ii)=[nodetemp(1:(groundnode-1),ii); 0; nodetemp((groundnode:end),ii)];
end   

%Solve forward problems on ALL measurement pairs
data=zeros(n_node+n_elec,totmeas);
%Put correct currents in the nodes in data matrix
count=0;
for k=1:n_stim
    %Find the measurement matrix for particular stimulation
    stimkmeasmatrix = stimstruc(k).meas_pattern; 
    %Loop over the measurement matrix
    for j=1:size(stimkmeasmatrix,1)
        %Find nodes to stimulate on measurement electrode
        stimkmeasjcur = stimkmeasmatrix(j,:);        
        %Add on to the counter
        count=count+1;
        %Put currents into RHS vector at the correct indexed node number
        for i=1:n_elec
            data(i+n_node,count)=stimkmeasjcur(i);
        end
    end
end
%Eliminate the ground node equation from data
data(groundnode,:)=[];
%Solve the sets of linear system (each column is solution of measure drive)
nodetemp=left_divide(At,data);
for ii=1:totmeas
    nodeunknownsmeas(:,ii)=[nodetemp(1:(groundnode-1),ii); 0; nodetemp((groundnode:end),ii)];
end   

%Calculate the structure of the COM of the normal and tangential space of 
%each electrode and the boundary normals.
elec_comp = calc_electrode_components(img.fwd_model);

%The Jacobian DOF's for each electrode is just 1 in 2D i.e. we find
%direction along vector (which is codim 1) and 2 in 3D i.e. direction along
%vector (which is codim 1).
dof_per_elec=nodedim-1; 
dofelec = length(elec_comp)*dof_per_elec;

%Initialise the Jacobian 
cntjac=0; J=zeros(totmeas,dofelec); 
for k=1:n_stim
    %Grab the stimulation nodal unknowns for this stimulation
    stimkunknowns = nodeunknownsfwd(:,k);

    %Loop over the measurement matrix
    for j=1:size(stimkmeasmatrix,1)
        %Add the measurement to the counter
        cntjac=cntjac+1;
 
        %Grab the measurement unknowns
        stimkmeasjunknowns=nodeunknownsmeas(:,cntjac);

        %Loop over each electrode and calculate the Jacobian
        %wr.t. average movement i.e. along each vector
        cntDOF=0; sens=zeros(1,dofelec);
        for i=1:length(elec_comp)
            for p=1:dof_per_elec   
                %Increase counter over electrode DOFs
                cntDOF=cntDOF+1;
                
                %Get the electrode measurement and stimulation voltage
                %(this only works for CEM so this is OK using n_node)
                UkjStim=stimkunknowns(n_node+i);
                UkjMeas=stimkmeasjunknowns(n_node+i);
            
                %Get the contact impedance
                zci = elecstruc(i).z_contact;            
                        
                %Get the end node numbers
                node_i_bound = elec_comp{i}.boundary_nodes; %2D should be two!!
            
                %Sample the nodes over electrode
                for ii=1:length(node_i_bound)
                    %Get the node number for this end point
                    node_n = node_i_bound(ii);
                    %Get the nodal measurement and stimulation voltage
                    ukjStim=stimkunknowns(node_n);
                    ukjMeas=stimkmeasjunknowns(node_n);    
        
                    %NOTE 3D NEED TO THINK MORE CAREFULLY i.e. WE NEED
                    %ANOTHER LOOP TO LOOK AT EACH COMPONENT OF MOVEMENT
                    
                    %Calculate (stim_elec - stim_node)*(meas_elec-meas_node)
                    sens(cntDOF) = sens(cntDOF) - ((UkjStim-ukjStim)*(UkjMeas-ukjMeas)/zci)* ...
                                       elec_comp{i}.boundary_nodes_dir(ii);                    
                end
            end
        end
        %Store the sensitivity in a matrix
        J(cntjac,:)=sens;  
    end
end

end
