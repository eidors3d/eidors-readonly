function J = jacobian_adjoint_higher_order(fwd_model,img)
%Find the Jacobian associated with an image (and forward model)
%Derivative of discretization method
% 
% Example (2D):
%  imdl = mk_common_model('c2C',16); img=mk_image(imdl.fwd_model,1);
%  img.fwd_model.solve = @fwd_solve_higher_order;
%  img.fwd_model.system_mat = @system_mat_higher_order;
%  img.fwd_model.jacobian = @jacobian_adjoint_higher_order;
%  
%  vve=[]; JJ4=[];
%  for i= 1:3; switch i;
%     case 1; img.fwd_model.approx_type = 'tri3'; % linear
%     case 2; img.fwd_model.approx_type = 'tri6'; % quadratic
%     case 3; img.fwd_model.approx_type = 'tri10'; % cubic;
%     end %switch
%     vv=fwd_solve(img);      vve(:,i)=vv.meas;
%     JJ=calc_jacobian(img);  JJ4(:,i)=JJ(4,:)';
%  end
%
% Example (3D):
%  imdl = mk_common_model('b3cr',16);  img=mk_image(imdl.fwd_model,1);
%  img.fwd_model.solve = @fwd_solve_higher_order;
%  img.fwd_model.system_mat = @system_mat_higher_order;
%  img.fwd_model.jacobian = @jacobian_adjoint_higher_order;
%  
%  vve=[]; JJ4=[];
%  for i= 1:2; switch i;
%     case 1; img.fwd_model.approx_type = 'tet4'; % linear
%     case 2; img.fwd_model.approx_type = 'tet10'; % quadratic
%     end %switch
%     vv=fwd_solve(img);      vve(:,i)=vv.meas;
%     JJ=calc_jacobian(img);  JJ4(:,i)=JJ(4,:)';
%  end

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_ADJOINT_HIGHER_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

if mdl_normalize(fwd_model)
     fwd_solve_data= fwd_solve( img );   
end

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
s_mat = calc_system_mat(img); At=s_mat.E; elemstiff=s_mat.elemstiff;
 
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

%Calculate Jacobian tensor - DE_{i,j,k} == dV_i,j / dS_k
%V_i,j - voltage change on electrode i for stim j
%S_k - conductivity change on element k
DE= zeros(nelecs,nstims,nelems);

%First step, we only want to pick off the ith electrode
zi2E(:,idx) = Node2Elec(:,idx)/At(idx,idx);

%SPEED UP HERE
%Factorise A = C'*S*C  - S diagonal conduc (C=system_mat_fields)
%We don't need extra multiplication in loop below
%only for piecewise linear FEM??
%
%zi2E= zeros(nelecs, nnodes+nelecs);
%zi2E(:,idx) = Node2Elec(:,idx)/At(idx,idx);
%zi2E=zi2E*FC'; sv=Fc*sv;

%Calculate the partial derivative matrix for kth change
for k=1:nelems    
    %kth element stiffness matrix, global nodes and index vector
    stiffk=elemstiff(k).elemstiff; nodesk=elemstruc(k,:); idx2=1:size(nodesk,2);
        
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
if mdl_normalize(fwd_model)
     data= fwd_solve_data; % must calculate first, because fwd_model is changed
     J= J ./ (data.meas(:)*ones(1,nelems));
end

%Negative Jacobian for injected currents??
J= -J;  
end

function do_unit_test
   tol = 1e-14;
   JJ = do_unit_test_2D(0); % not normalized
   JJ_ref = -1e-4*[
   2.440415272063380   2.705754096983918   2.721135010947898
   2.578623854199123   2.327923064064513   2.342086727271722
   1.438743206711758   1.333580385504260   1.337599904647092
   1.300534624576032   1.650702059478922   1.659896278693538];

   unit_test_cmp('2D: 1st order',JJ(1:4,1),JJ_ref(1:4,1),tol);
   unit_test_cmp('2D: 2nd order',JJ(1:4,2),JJ_ref(1:4,2),tol);
   unit_test_cmp('2D: 3rd order',JJ(1:4,3),JJ_ref(1:4,3),tol);

   [JJ1,vve]= do_unit_test_2D(1); for i=1:3
   unit_test_cmp('2D: (normalize)',JJ1(:,i),JJ(:,i)/vve(4,i),tol);
   end

   JJ = do_unit_test_3D(0);
   JJ_ref = -1e-5*[
   1.246064580179371   1.585061706092707
   1.332632578853691   1.354929239220232
   0.712061825721561   0.443297935900921
   0.625493827047241   0.604174950085724];
   [JJ1,vve]= do_unit_test_3D(1); for i=1:2
   unit_test_cmp('3D: (normalize)',JJ1(:,i),JJ(:,i)/vve(4,i),tol);
   end

   unit_test_cmp('3D: 1st order',JJ(1:4,1),JJ_ref(1:4,1),tol);
   unit_test_cmp('3D: 2nd order',JJ(1:4,2),JJ_ref(1:4,2),tol);

end
function [JJ4,vve]=do_unit_test_2D(normalize_flag)
   imdl = mk_common_model('c2C',16); img = mk_image(imdl.fwd_model,1);
   img.fwd_model.normalize_measurements = normalize_flag;
   vv=fwd_solve(img);      v0e=vv.meas;
   JJ=calc_jacobian(img);  J04=JJ(4,:)';

   %High-order EIDORS solver %Change default eidors solvers
   img.fwd_model.solve = @fwd_solve_higher_order;
   img.fwd_model.system_mat = @system_mat_higher_order;
   img.fwd_model.jacobian = @jacobian_adjoint_higher_order;

   vve=[]; JJ4=[];
   for i= 1:3; switch i;
      case 1; img.fwd_model.approx_type = 'tri3'; % linear
      case 2; img.fwd_model.approx_type = 'tri6'; % quadratic
      case 3; img.fwd_model.approx_type = 'tri10'; % cubic;
      end %switch
      vv=fwd_solve(img);      vve(:,i)=vv.meas;
      JJ=calc_jacobian(img);  JJ4(:,i)=JJ(4,:)';
   end

   subplot(321);
   plot([v0e,vve,(v0e*[1,1,1]-vve)*10]);
   legend('Default','linear','quadratic','cubic','(1-0)x10','(2-0)x10','(3-0)x10');
   xlim([1,100]);

   imgJJ=img; imgJJ.elem_data = JJ4;
   imgJJ.show_slices.img_cols = 3;

   subplot(323); show_slices(imgJJ); eidors_colourbar(imgJJ);

   imgJJ.elem_data = JJ4 - J04*[1,1,1];
   subplot(325); show_slices(imgJJ); eidors_colourbar(imgJJ);
end
function [JJ4,vve]=do_unit_test_3D(normalize_flag)
   imdl = mk_common_model('b3cr',16); img = mk_image(imdl.fwd_model,1);
   img.fwd_model.normalize_measurements = normalize_flag;
   vv=fwd_solve(img);      v0e=vv.meas;
   JJ=calc_jacobian(img);  J04=JJ(4,:)';

   %High-order EIDORS solver %Change default eidors solvers
   img.fwd_model.solve = @fwd_solve_higher_order;
   img.fwd_model.system_mat = @system_mat_higher_order;
   img.fwd_model.jacobian = @jacobian_adjoint_higher_order;

   vve=[]; JJ4=[];
   for i= 1:2; switch i;
      case 1; img.fwd_model.approx_type = 'tet4'; % linear
      case 2; img.fwd_model.approx_type = 'tet10'; % quadratic
      end %switch
      vv=fwd_solve(img);      vve(:,i)=vv.meas;
      JJ=calc_jacobian(img);  JJ4(:,i)=JJ(4,:)';
   end

   subplot(322);
   plot([v0e,vve,(v0e*[1,1]-vve)*10]);
   legend('Default','linear','quadratic','(1-0)x10','(2-0)x10');
   xlim([1,100]);

   imgJJ=img; imgJJ.elem_data = JJ4;
   imgJJ.show_slices.img_cols = 2;

   level = [inf,inf,0.3];
   subplot(324); show_slices(imgJJ,level); eidors_colourbar(imgJJ);

   imgJJ.elem_data = JJ4 - J04*[1,1];
   subplot(326); show_slices(imgJJ,level); eidors_colourbar(imgJJ);
end
