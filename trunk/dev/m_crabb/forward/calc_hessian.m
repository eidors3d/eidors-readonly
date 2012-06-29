function H = calc_hessian(fwd_model,img)
%Find the Hessian associated with an image (and forward model)
%Second derivative of discretization method
%This is just a cached version of code - mc_calc_hessian

cache_obj= {fwd_model,img};
H = eidors_obj('get-cache', cache_obj, 'hessian');
if ~isempty(H)
            eidors_msg('hessian: using cached value', 3);        
        return
end

[H]= mc_calc_hessian(fwd_model,img);

eidors_obj('set-cache', cache_obj, 'hessian',H);
eidors_msg('hessian: setting cached value', 3);

end

function H = mc_calc_hessian(fwd_model,img)
%Find the Hessian associated with an image (and forward model)
%Second derivative of discretization method

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end

%Modify the forward model to be of my type
%%fwd_model = mc_fem_modify(fwd_model); img.fwd_model=fwd_model;
[bound,elem,nodes] = fem_modify(fwd_model); 
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

%Calculate Hessian tensor - D2E_{i,j,k,l} == d2V_i,j / dS_k dS_l
%V_i,j - voltage change on electrode i for stim j
%S_k/S_l - conductivity change on element k and elemen l
D2E= zeros(nelecs,nstims,nelems,nelems);

%First step, we only want to pick off the ith electrode
zi2E(:,idx) = Node2Elec(:,idx)/At(idx,idx);

%Calculate the partial derivative matrix for kth change
for k=1:nelems
        %Get the kth element global nodes and stiffness
        stiffk=elemstiff(k).stiff; nodesk=elem(k).nodes; idxk=1:size(nodesk,2);
        
        %Create the FEM derivative matrix and multiply by inverse
        dA_dSk=dA_zero; dA_dSk(nodesk(idxk),nodesk(idxk))=stiffk(idxk,idxk);
        dA_dSk2E=dA_dSk(:,idx)/At(idx,idx);
        
    for l=1:nelems
        %Get the lth element global nodes and stiffness
        stiffl=elemstiff(l).stiff; nodesl=elem(l).nodes; idxl=1:size(nodesl,2);
        
        %Create the FEM derivative matrix
        dA_dSl=dA_zero; dA_dSl(nodesk(idxl),nodesk(idxl))=stiffl(idxl,idxl);
    
        %Now form product with solution
        D2E(:,:,k,l) = zi2E(:,idx)*dA_dSk2E(idx,:)*dA_dSl(idx,idx)*nodeunknownsfwd(idx,:);
    end
end

%Calculate Hessian tensor (measurement patterns specified here)
cnthes=0; H=zeros(nmeass,nelems,nelems);
for j=1:nstims   
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   D2Ej = reshape( D2E(:,j,:,:), nelecs, nelems, nelems);
   %Multiplication in loop, matlab doesnt like tensor multiplication
   for kk=1:nelems %Loop over elements
       D2Ej_kk=reshape(D2Ej(:,kk,:),nelecs,nelems);
       H(cnthes+(1:n_meas),kk,:) = meas_pat*D2Ej_kk;
   end
   cnthes = cnthes + n_meas;
end 

%Get the Jacobian and normalize measurements (if field exists)
if isfield(fwd_model,'normalize_measurements')
%    data=mc_fwd_solve( img );   
%    J= J ./ (data.meas(:)*ones(1,nelems));
end

end
