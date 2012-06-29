function J = jacobian_adjoint_higher_order(fwd_model,img)
%Find the Jacobian associated with an image (and forward model)
%Derivative of discretization method

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Modify the forward model to be of my type
%%fwd_model = mc_fem_modify(fwd_model); img.fwd_model=fwd_model;
[bound,elem,nodes] = fem_1st_to_higher_order(fwd_model); 
fwd_model.bound=bound; fwd_model.elem=elem; fwd_model.nodes=nodes;
img.fwd_model=fwd_model;

%Calculate the total stiffness matrix and elemental stiffness matrices
s_mat = calc_system_mat(fwd_model,img); At=s_mat.E;

%Get the stiffness matrix structure
elemstiff = mc_calc_elem_stiff(fwd_model); %CHANGE THIS!!!!!!!!!!!!!
 
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

%Complete or Point? - Check first electrode and change index vector of 
%'node' number corresponding to electrode
elecnode=zeros(1,nelecs);
if(size(elecstruc(1).nodes,2)==1 && size(elecstruc(1).nodes,1)==1) %POINT
    %Initialise node to electrode matrix
    Node2Elec=sparse(nelecs,nnodes);
    for i=1:nelecs
        %Assign electrode index at correct node
        elecnode(i)=elecstruc(i).nodes;
        Node2Elec(i,elecnode(i))=1;
    end
    %Assign a matrix for derivative of FEM w.r.t conduc
    dA_zero=sparse(nnodes,nnodes);
    
    %Assign correct size unknowns and right hand side matrix (forward)
    datafwd=zeros(nnodes,nstims); 
    nodeunknownsfwd=zeros(nnodes,nstims); 
else
    %Initialise node to electrode matrix
    Node2Elec=sparse(nelecs,nnodes+nelecs);
    for i=1:nelecs
        %Assign electrode at bottom of list
        elecnode(i)=nnodes+i;
        Node2Elec(i,elecnode(i))=1;
    end
    
    %Assign a matrix for derivative of FEM w.r.t conduc
    dA_zero=sparse(nnodes+nelecs,nnodes+nelecs);
        
    %Assign correct size unknowns and right hand side matrix (forward)
    datafwd=zeros(nnodes+nelecs,nstims); 
    nodeunknownsfwd=zeros(nnodes+nelecs,nstims); 
end

%Loop over stimulations and assign current matrix
%CHANGE THIS BY USING NODE2ELEC MATRIX!!!!
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

%Calculate Jacobian tensor - DE_{i,j,k} == dV_i,j / dS_k
%V_i,j - voltage change on electrode i for stim j
%S_k - conductivity change on element k
DE= zeros(nelecs,nstims,nelems);

%First step, we only want to pick off the ith electrode
zi2E(:,idx) = Node2Elec(:,idx)/At(idx,idx);

%Calculate the partial derivative matrix for kth change
for k=1:nelems    
    %Get the kth element stiffness matrix and global nodes
    stiffk=elemstiff(k).stiff; nodesk=elem(k).nodes;
    
    %Get an index vector for the 
    idx2=1:size(nodesk,2);
        
    %Create the FEM derivative matrix
    dA_dSk=dA_zero; dA_dSk(nodesk(idx2),nodesk(idx2))=stiffk(idx2,idx2);

    %Now form product with solution
    DE(:,:,k) = zi2E(:,idx)*dA_dSk(idx,idx)*nodeunknownsfwd(idx,:);
end

%Calculate Jacobian matrix (measurement patterns specified here)
cntjac=0; J=zeros(nmeass,nelems);
for j=1:nstims   
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   DEj = reshape( DE(:,j,:), nelecs, nelems);
   J( cntjac+(1:n_meas),: ) = meas_pat*DEj;
   cntjac = cntjac + n_meas;
end; 

%Get the Jacobian and normalize measurements (if field exists)
if isfield(fwd_model,'normalize_measurements')
%    data=mc_fwd_solve( img );   
%    J= J ./ (data.meas(:)*ones(1,nelems));
end

%Negative Jacobian for injected currents??
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
[weight,xcoord,ycoord,zcoord]=element_gauss_points(eletype);

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
            (jacobianelem\element_d_shape_function(eletype,xcoord(kk),ycoord(kk),zcoord(kk)))'* ...
            (jacobianelem\element_d_shape_function(eletype,xcoord(kk),ycoord(kk),zcoord(kk)))*magjacelem;
    end
    %This is element stiffness matrix (multiply by conductivity????)
    elemstiff(i).stiff=Ammat; %elemstiff(i).stiff=Ammat*img.elem_data(i); 
end

end
