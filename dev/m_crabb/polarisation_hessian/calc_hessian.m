function H = calc_hessian(fwd_model,img,elem_list)
%Find the Hessian associated with an image (and forward model)
%Second derivative of discretization method d2F_dskdsl
%3rd input elem_list - indices of elementn we want Hessian for
%M Crabb - 29.06.2012
%
%TODO - (i) Use transpose to speed up calc i.e. D2H = H+H'?????

cache_obj= {fwd_model,img};
H = eidors_obj('get-cache', cache_obj, 'hessian');
if ~isempty(H)
            eidors_msg('hessian: using cached value', 3);        
        return
end


[H]= mc_calc_hessian(fwd_model,img,elem_list);

eidors_obj('set-cache', cache_obj, 'hessian',H);
eidors_msg('hessian: setting cached value', 3);

end

function H = mc_calc_hessian(fwd_model,img,elem_list)
%Find the Hessian associated with an image (and forward model)
%Second derivative of discretization method

%If function called only with image, extract forward model
if(nargin==1)
    img=fwd_model; fwd_model=img.fwd_model;
end


%Modify the forward model to be of my type
%%fwd_model = mc_fem_modify(fwd_model); img.fwd_model=fwd_model;
%[bound,elem,nodes] = fem_1st_to_higher_order(fwd_model); 
%fwd_model.bound=boundary; fwd_model.elem=elems; fwd_model.nodes=nodes;
%img.fwd_model=fwd_model;

%Calculate the total stiffness matrix and elemental stiffness matrices
s_mat = system_mat_higher_order(fwd_model,img); 
At=s_mat.E; elemstiff=s_mat.elemstiff;
 
%Find electrode stucture and no.of electrodes 
%Find stim strucutre and no. stimulations
%Find node structure and find no.nodes 
%Find element structure and create vector of length no. elements
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 
elemstruc=fwd_model.elems; nelems=size(elemstruc,1); 

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
D2E= zeros(nelecs,nstims,length(elem_list),length(elem_list));
%First step, we only want to pick off the ith electrode
zi2E(:,idx) = Node2Elec(:,idx)/At(idx,idx);

   
%Calculate the partial derivative matrix for kth change
for idxkk=1:length(elem_list)
        %Select kth element
        k=elem_list(idxkk);
    
        %Get the kth element global nodes and stiffness
        stiffk=elemstiff(k).elemstiff; nodesk=elemstruc(k,:); idxk=1:size(nodesk,2);
        
        %Create the FEM derivative matrix and multiply by inverse
        dA_dSk=dA_zero; dA_dSk(nodesk(idxk),nodesk(idxk))=stiffk(idxk,idxk);
        dA_dSk2E=dA_dSk(:,idx)/At(idx,idx);
        
    parfor idxll=1:length(elem_list)
        %Select kth element
        l=elem_list(idxll);      
        
        %Get the lth element global nodes and stiffness
        stiffl=elemstiff(l).elemstiff; nodesl=elemstruc(l,:); idxl=1:size(nodesl,2);
        
        %Create the FEM derivative matrix
        dA_dSl=dA_zero; dA_dSl(nodesl(idxl),nodesl(idxl))=stiffl(idxl,idxl);
        dA_dSl2E=dA_dSl(:,idx)/At(idx,idx);

        %Now form product with solution
        %d2Udldk = Z *( A^{-1}DA(k)A^{-1}DA(l)A^{-1} + A^{-1}DA(l)A^{-1}DA(k)A^{-1} )f =    
        D2E(:,:,idxkk,idxll) = zi2E(:,idx)*dA_dSk2E(idx,:)*dA_dSl(idx,idx)*nodeunknownsfwd(idx,:) + ...
                               zi2E(:,idx)*dA_dSl2E(idx,:)*dA_dSk(idx,idx)*nodeunknownsfwd(idx,:);
    end
    fprintf(1,'Hessian completed for element %i of %i\n',idxkk,length(elem_list));
end

%Calculate Hessian tensor (measurement patterns specified here)
cnthes=0; H=zeros(nmeass,length(elem_list),length(elem_list));
for j=1:nstims   
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   D2Ej = reshape( D2E(:,j,:,:), nelecs, length(elem_list), length(elem_list) );
   %Multiplication in loop, matlab doesnsize(D2Et like tensor multiplication
   for idxkk=1:length(elem_list) %Loop over elements
       D2Ej_kk=reshape(D2Ej(:,idxkk,:),nelecs,length(elem_list));
       H(cnthes+(1:n_meas),idxkk,:) = meas_pat*D2Ej_kk;
   end
   cnthes = cnthes + n_meas;
end 

end