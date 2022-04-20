function[data] = fwd_solve_higher_order(fwd_model,img)
%Solve for voltages (nodes/electrodes) for a forward model. 
%M Crabb - 29.06.2012
% 
% Example (2D):
%  imdl = mk_common_model('c2C',16); img=mk_image(imdl.fwd_model,1);
%  img.fwd_model.solve = @fwd_solve_higher_order;
%  img.fwd_model.system_mat = @system_mat_higher_order;
%  
%  vve=[]; JJ4=[];
%  for i= 1:3; switch i;
%     case 1; img.fwd_model.approx_type = 'tri3'; % linear
%     case 2; img.fwd_model.approx_type = 'tri6'; % quadratic
%     case 3; img.fwd_model.approx_type = 'tri10'; % cubic;
%     end %switch
%     vv=fwd_solve(img);      vve(:,i)=vv.meas;
%  end
%
% Example (3D):
%  imdl = mk_common_model('b3cr',16);  img=mk_image(imdl.fwd_model,1);
%  img.fwd_model.solve = @fwd_solve_higher_order;
%  img.fwd_model.system_mat = @system_mat_higher_order;
%  
%  vve=[]; JJ4=[];
%  for i= 1:2; switch i;
%     case 1; img.fwd_model.approx_type = 'tet4'; % linear
%     case 2; img.fwd_model.approx_type = 'tet10'; % quadratic
%     end %switch
%     vv=fwd_solve(img);      vve(:,i)=vv.meas;
%  end

if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling FWD_SOLVE_HIGHER_ORDER with two arguments is deprecated and will cause' ...
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
    %We need to update fwd_model of img too for system_mat
    img.fwd_model=fwd_model;
end

%Calculate the total stiffness matrix
s_mat = calc_system_mat(img); At=s_mat.E;

%Find electrode stucture and no.of electrodes and initialize vector
%Find stim strucutre and no. stimulations
%Find node structure and find no.nodes 
elecstruc=fwd_model.electrode; nelecs=size(elecstruc,2);
stimstruc=fwd_model.stimulation; nstims=size(stimstruc,2); 
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
nodeunknowns(idx,:)=left_divide(At(idx,idx),rhsdata(idx,:));


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
    % Use conj on meas_pat: see Chap 5 Adler&Holder 2021
    meas_pat=conj(stimstruc(ii).meas_pattern); %Measurement patterns
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


end

function do_unit_test
   tol = 1e-13;
   vve = do_unit_test_2D;
   vve_ref = [
   0.898225115241117   0.921510761628436   0.928516596253320
   0.406398239528486   0.412406966923347   0.413938179536629
   0.248067950415984   0.250212111398160   0.250609888163540
   0.179592593985273   0.179643951982781   0.179734695991927];

   unit_test_cmp('2D: 1st order',vve(1:4,1),vve_ref(1:4,1),tol);
   unit_test_cmp('2D: 2nd order',vve(1:4,2),vve_ref(1:4,2),tol);
   unit_test_cmp('2D: 3rd order',vve(1:4,3),vve_ref(1:4,3),tol);
   vve = do_unit_test_3D;
   vve_ref = [
   1.404189968566952   1.410900674463290
   0.403207625837809   0.402992578774667
   0.198517193844915   0.201211971071238
   0.133852904079284   0.133841105904217];
   unit_test_cmp('3D: 1st order',vve(1:4,1),vve_ref(1:4,1),tol);
   unit_test_cmp('3D: 2nd order',vve(1:4,2),vve_ref(1:4,2),tol);

end
function [vve] = do_unit_test_2D
   imdl = mk_common_model('c2C',16); img = mk_image(imdl.fwd_model,1);
   vv=fwd_solve(img);      v0e=vv.meas;

   %High-order EIDORS solver %Change default eidors solvers
   img.fwd_model.solve = @fwd_solve_higher_order;
   img.fwd_model.system_mat = @system_mat_higher_order;

   vve=[]; JJ4=[];
   for i= 1:3; switch i;
      case 1; img.fwd_model.approx_type = 'tri3'; % linear
      case 2; img.fwd_model.approx_type = 'tri6'; % quadratic
      case 3; img.fwd_model.approx_type = 'tri10'; % cubic;
      end %switch
      vv=fwd_solve(img);      vve(:,i)=vv.meas;
   end

   subplot(321);
   plot([v0e,vve,(v0e*[1,1,1]-vve)*10]);
   legend('Default','linear','quadratic','cubic','(1-0)x10','(2-0)x10','(3-0)x10');
   xlim([1,100]);
end
function vve=do_unit_test_3D;
   imdl = mk_common_model('b3cr',16); img = mk_image(imdl.fwd_model,1);
   vv=fwd_solve(img);      v0e=vv.meas;

   %High-order EIDORS solver %Change default eidors solvers
   img.fwd_model.solve = @fwd_solve_higher_order;
   img.fwd_model.system_mat = @system_mat_higher_order;

   vve=[]; JJ4=[];
   for i= 1:2; switch i;
      case 1; img.fwd_model.approx_type = 'tet4'; % linear
      case 2; img.fwd_model.approx_type = 'tet10'; % quadratic
      end %switch
      vv=fwd_solve(img);      vve(:,i)=vv.meas;
   end

   subplot(322);
   plot([v0e,vve,(v0e*[1,1]-vve)*10]);
   legend('Default','linear','quadratic','(1-0)x10','(2-0)x10');
   xlim([1,100]);

end
