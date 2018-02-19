function data = calc_neumann_func_nodal(fwd_model,img)

% \nabla \cdot (A \nabla N (x,y) )= -\delta_{y}
% \nu \cdot (A \nabla N (x,y))|_{\partial \Omega} = -1/|\partial \Omega|
% \int_{\partial \Omega} \nu \cdot (A \nabla N (x,y))|_{\partial \Omega}
%
% Ammari/Kang - Review - pg.36 

% NB1 - I think this is independent of CEM i.e. this comes in from the
% other simulated data???
%
% N is size nodes x elements i.e the columns are the Greens function (which
% are essentially voltage potentials, and hence node based) associated with
% a delta source for each finite element

%Data - Neumann Greens function. For FE system this is dim n_node*n_node

%% JUST LINEAR ELEMENT ATM

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling FWD_SOLVE_HIGHER_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end

fwd_model= img.fwd_model;
n_dim = size(fwd_model.nodes,2);

%Modify the forward model to be of my type
if ~isfield(fwd_model,'approx_type')    || ...
   strcmp(fwd_model.approx_type,'tri3') || ...
   strcmp(fwd_model.approx_type,'tet4')   
    %Do nothing
else
    [bound,elem,nodes] = fem_1st_to_higher_order(fwd_model);
    fwd_model.boundary=bound; fwd_model.elems=elem; fwd_model.nodes=nodes;
    %We need to update fwd_model of img too for system_mat
    img.fwd_model=fwd_model;
end


%Find electrode stucture and no.of electrodes and initialize vector
%Find stim strucutre and no. stimulations
%Find node structure and find no.nodes 
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
elemstruc=fwd_model.elems; nelems=size(elemstruc,1);
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 
nodestruc=fwd_model.nodes; nnodes=size(nodestruc,1); 
boundstruc=img.fwd_model.boundary;  nbounds = size(boundstruc,1);

%Calculate the total stiffness matrix
s_mat = stiffness_mat_higher_order(img); At=s_mat.E(1:nnodes,1:nnodes);

%Just using pcw linear
eletype = 'tri3';

%Gauss quadrature points to perform numerical integration
[weight,xcoord,ycoord]=gauss_points(1,6);
for k=size(weight,2):-1:1
    phi(:,k) = boundary_shape_function(eletype,xcoord(k),ycoord(k))';
end
nodestruc=img.fwd_model.nodes; 

%Now initialise the rhs_bound vector
rhs_bound=zeros(nnodes,1);

%%Compute the length/area of boundary
%bound_volume = 0;
%for i=1:n_bound;
%    thisb=nodestruc(boundstruc(i,:),:);    
%
%    if(n_dim==2)    
%        diff21=thisb(2,:)-thisb(1,:);        
%        bound_volume_i = norm(diff21);
%
%    end
[bdy_idx,bdy_area] = find_electrode_bdy(boundstruc(:,1:n_dim),nodestruc,unique(boundstruc));
bound_area = sum(bdy_area);

nodedim=size(nodestruc,2);

%An extra row equation for constrain that integral is 0
Aw = zeros(1,nnodes);

%Integrate RHS vector over boundary
for i=1:nbounds
    %Jacobian of mapping
    boundnodes=boundstruc(i,:);    
    thisb=nodestruc(boundnodes,:);    

    %Find the magnitude Jacobian of the mapping in 2D/3D
    %NB:Scalings are consistent with reference element shape
    if(nodedim==2)
        %Jacobian = 0.5*|(x2-x1)| (x1,x2 vector of coords)
        diff21=thisb(2,:)-thisb(1,:);
        magjacbound=0.5*sqrt(diff21(1)^2+diff21(2)^2);
    elseif(nodedim==3)
        %Jacobian = |(x3-x1)x(x3-x2)| (x1,x2,x3 vector of coords)
        diffprod=cross(thisb(3,:)-thisb(1,:),thisb(3,:)-thisb(2,:));
        magjacbound=sqrt(diffprod(1)^2+diffprod(2)^2+diffprod(3)^2);
    end
    %Initilaise rhs vector of boundary and integral 0 condition 
    fmat1=0;       Awmat=0;
    for k=1:size(weight,2)
        %Find global coordinates of wiehgt point        
        fmat1 = fmat1 + weight(k)*(phi(:,k))'*(-1/bound_area)*magjacbound;    
        Awmat = Awmat + weight(k)*(phi(:,k))'*magjacbound;           
    end

    %Assemble the matrices
    Aw(1,boundnodes,1) = Aw(1,boundnodes) + Awmat;                
    rhs_bound(boundnodes)=rhs_bound(boundnodes)+ fmat1';
    
end

rhs_total = zeros(nnodes+1,nnodes);
%rhs_total(1,:)
for i=1:nnodes %For each triangle
   rhs_total(1:nnodes,i) = rhs_bound;
   
   %This is 1 at the ith node
   rhs_total(i,i) = rhs_total(i,i) +1;    
end
rhs_total(nnodes+1,:)=0;

%Add extra equation i.e. At-> (n+1) x n dimension full column rank
At=[At;Aw];

%Solve the simulated linear system with index
nodeunknowns=left_divide(At,rhs_total);

data = nodeunknowns;








%{
%Find electrode voltages and store in matrix
%Calculate electrode voltages using index vector elecnode(i)
velec=zeros(nelecs,nstims);
for i=1:nelecs
    %This is the indexed entries in nodeunknowns
    velec(i,:)=nodeunknowns(elecnode(i),:);
end

%Get the measured voltages 
vmeaselec=zeros(nmeass,1); idx=0;
for ii=1:nstims
    meas_pat=stimstruc(ii).meas_pattern; %Measurement patterns
    n_meas=size(meas_pat,1); %Number of measures
    vmeaselec(idx + (1:n_meas) ) = meas_pat*velec(:,ii); %Diff data
    idx=idx+n_meas; %Increase counter
end

%Return the electrode voltages in data structure
data.meas= vmeaselec;
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_higher_order';
data.quantity = 'voltage';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = nodeunknowns(1:nnodes,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = nodeunknowns;             % all, including CEM nodes
end; end
%}
