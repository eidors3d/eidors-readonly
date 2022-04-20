function data =fwd_solve_1st_order(fwd_model, img)
% FWD_SOLVE_1ST_ORDER: data= fwd_solve_1st_order( img)
% First order FEM forward solver
% Input:
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)
%    img.fwd_solve.get_elec_curr= 1 (data.elec_curr = current on electrodes)
%
% Model Reduction: use precomputed fields to reduce the size of
%    the forward solution. Nodes which are 1) not used in the output
%    (i.e. not electrodes) 2) all connected to the same conductivity via
%    the c2f mapping are applicable.
% see: Model Reduction for FEM Forward Solutions, Adler & Lionheart, EIT2016
%
%    img.fwd_model.model_reduction = @calc_model_reduction;
%       where the functionputs a struct with fields: main_region, regions
%       OR
%    img.fwd_model.model_reduction.main_region = vector, and 
%    img.fwd_model.model_reduction.regions = struct

% (C) 1995-2017 Andy Adler. License: GPL version 2 or version 3
% $Id$

% correct input paralemeters if function was called with only img
if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling FWD_SOLVE_1ST_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

img = data_mapper(img);
if ~ismember(img.current_params, supported_params)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'FWD_SOLVE_1ST_ORDER',img.current_params);
end
% all calcs use conductivity
img = convert_img_units(img, 'conductivity');

pp= fwd_model_parameters( fwd_model, 'skip_VOLUME' );
pp = set_gnd_node(fwd_model, pp);
s_mat= calc_system_mat( img );

if isfield(fwd_model,'model_reduction')
   [s_mat.E, main_idx, pp] = mdl_reduction(s_mat.E, ...
           img.fwd_model.model_reduction, img, pp );
else
   pp.mr_mapper = 1:size(s_mat.E,1);
end

% Normally EIT uses current stimulation. In this case there is
%  only a ground node, and this is the only dirichlet_nodes value.
%  In that case length(dirichlet_nodes) is 1 and the loop runs once
% If voltage stimulation is done, then we need to loop on the
%  matrix to calculate faster. 
[dirichlet_nodes, dirichlet_values, neumann_nodes, has_gnd_node]= ...
         find_dirichlet_nodes( fwd_model, pp );

v= full(horzcat(dirichlet_values{:})); % Pre fill in matrix
for i=1:length(dirichlet_nodes)
   idx= 1:size(s_mat.E,1);
   idx( dirichlet_nodes{i} ) = [];
   % If all dirichlet patterns are the same, then calc in one go
   if length(dirichlet_nodes) == 1; rhs = 1:size(pp.QQ,2);
   else                           ; rhs = i; end
   v(idx,rhs)= left_divide( s_mat.E(idx,idx), ...
                     neumann_nodes{i}(idx,:) - ...
                     s_mat.E(idx,:)*dirichlet_values{i},fwd_model);
end

% If model has a ground node, check if current flowing in this node
if has_gnd_node
   Ignd = s_mat.E(dirichlet_nodes{1},:)*v;
   Irel = Ignd./sum(abs(pp.QQ)); % relative 
   if max(abs(Irel))>1e-6
      eidors_msg('%4.5f%% of current is flowing through ground node. Check stimulation pattern', max(abs(Irel))*100,1);
   end
end

% calc voltage on electrodes

% This is horribly inefficient, override
% v_els= pp.N2E * v;
idx = find(any(pp.N2E));
v_els= pp.N2E(:,idx) * v(idx,:);


% create a data structure to return
data.meas= meas_from_v_els(v_els, pp);
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_1st_order';
try; if img.fwd_solve.get_all_meas == 1
   outmap = pp.mr_mapper(1:pp.n_node);
   data.volt(outmap,:) = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt(pp.mr_mapper,:) = v;                % all, including CEM nodes
end; end
try; if img.fwd_solve.get_elec_curr== 1
%  data.elec_curr = pp.N2E * s_mat.E * v;
   idx = find(any(pp.N2E));
   data.elec_curr = pp.N2E(:,idx) * s_mat.E(idx,:) * v;
end; end


% has_gnd_node = flag if the model has a gnd_node => can warn if current flows
function [dirichlet_nodes, dirichlet_values, neumann_nodes, has_gnd_node]= ...
            find_dirichlet_nodes( fwd_model, pp );
   fnanQQ = isnan(pp.QQ);
   lstims = size(pp.QQ,2);
   % Can't use any(...) because if does implicit all
   if any(any(fnanQQ))
      has_gnd_node = 0; % no ground node is specified
      % Are all dirichlet_nodes the same

      % Don't use all on sparse, it will make them full
      % Check if all rows are the same
      if ~any(any(fnanQQ(:,1)*ones(1,lstims) - fnanQQ,2))
         dirichlet_nodes{1} = find(fnanQQ(:,1));
         dirichlet_values{1} = sparse(size(pp.N2E,2), size(fnanQQ,2));
         dirichlet_values{1}(fnanQQ) = pp.VV(fnanQQ);
         neumann_nodes{1} = pp.QQ;
         neumann_nodes{1}(fnanQQ) = 0;
      else % one at a time
         for i=1:size(fnanQQ,2)
            fnanQQi= fnanQQ(:,i);
            if any(fnanQQi)
               dirichlet_nodes{i} = find(fnanQQi);
               dirichlet_values{i} = sparse(size(pp.N2E,2), 1);
               dirichlet_values{i}(fnanQQi) = pp.VV(fnanQQi,i);
               neumann_nodes{i} = pp.QQ(:,i);
               neumann_nodes{i}(fnanQQi) = 0;
            elseif isfield(pp,'gnd_node')
               dirichlet_nodes{i} = pp.gnd_node;
               dirichlet_values{i} = sparse(size(pp.N2E,2), 1);
               neumann_nodes{1}   = pp.QQ(:,i);
               has_gnd_node= 1;
            else
               error('no required ground node on model');
            end
         end
      end
   elseif isfield(pp,'gnd_node')
      dirichlet_nodes{1} = pp.gnd_node;
      dirichlet_values{1} = sparse(size(pp.N2E,2), size(fnanQQ,2));
      neumann_nodes{1}   = pp.QQ;
      has_gnd_node= 1;
   else
      error('no required ground node on model');
   end

function pp = set_gnd_node(fwd_model, pp);
   if isfield(fwd_model,'gnd_node');
      pp.gnd_node = fwd_model.gnd_node;
   else
      % try to find one in the model center
      ctr =  mean(fwd_model.nodes,1);
      d2  =  sum(bsxfun(@minus,fwd_model.nodes,ctr).^2,2);
      [~,pp.gnd_node] = min(d2);
      eidors_msg('Warning: no ground node found: choosing node %d',pp.gnd_node(1),1);
   end

function vv = meas_from_v_els( v_els, pp)
   try
% Was 1.82s
%        % measured voltages from v
%    %   vv = zeros( pp.n_meas, 1 );
%        idx=0;
%        for i=1:length(stim)
%           meas_pat= stim(i).meas_pattern;
%           n_meas  = size(meas_pat,1);
%           vv( idx+(1:n_meas) ) = meas_pat*v_els(:,i);
%           idx= idx+ n_meas;
%        end

% This code replaced the previous - Nov 4, 2013
% Now 0.437s
% Why is it faster??

       [n_elec,n_stim] = size(v_els);

       vv = pp.v2meas * v_els(:);
   catch err
      if strcmp(err.identifier, 'MATLAB:innerdim');
          error(['measurement pattern not compatible with number' ...
                  'of electrodes for stimulation patter %d'],i);
      else
          rethrow(err);
      end
   end

   
function v2meas = get_v2meas(n_elec,n_stim,stim)
    v2meas = sparse(n_elec*n_stim,0);
    for i=1:n_stim
        meas_pat= stim(i).meas_pattern;
        n_meas  = size(meas_pat,1);
        v2meas((i-1)*n_elec + 1: i*n_elec,end+(1:n_meas)) = meas_pat.';
    end


function [E, m_idx, pp] = mdl_reduction(E, mr, img, pp);
   % if mr is a string we assume it's a function name
   if isa(mr,'function_handle') || ischar(mr)
      mr = feval(mr,img.fwd_model);
   end
   % mr is now a struct with fields: main_region, regions
   m_idx = mr.main_region;
   E = E(m_idx, m_idx);
   for i=1:length(mr.regions)
      invEi=   mr.regions(i).invE;
% FIXME:!!! data_mapper has done the c2f. But we don't want that here.
%  kludge is to reach into the fine model field. This is only ok because
%  model_reduction is only valid if one parameter describes each field
      field = mr.regions(i).field; 
      field = find(img.fwd_model.coarse2fine(:,field));
      field = field(1); % they're all the same - by def of model_reduction
      sigma = img.elem_data(field);
      E = E - sigma*invEi;
   end

   % Adjust the applied current and measurement matrices
   pp.QQ = pp.QQ(m_idx,:);
   pp.VV = pp.VV(m_idx,:);
   pp.N2E= pp.N2E(:,m_idx);
   pp.mr_mapper = cumsum(m_idx); %must be logical
   pp.gnd_node = pp.mr_mapper(pp.gnd_node);
   if pp.gnd_node==0
      error('model_reduction removes ground node');
   end
   
        
function unit_test_voltage_stims;
   stim = zeros(16,1); volt=stim; stim([1,4]) = NaN; volt([1,4]) = [1,2];
   stimulv.stim_pattern = stim;
   stimulv.volt_pattern = volt;
   stimulv.meas_pattern = [1,0,0,-1,zeros(1,12)];

   img = mk_image( mk_common_model('a2c2',16),1);
   img.fwd_model.stimulation = stimulv;
   img.fwd_solve.get_all_nodes = 1;
   vh = fwd_solve_1st_order(img);
   unit_test_cmp('a2c2 Vstim #1', vh.meas, diff(volt(1:2,1)), 1e-14);
   unit_test_cmp('a2c2 Vstim #2', vh.volt(27+[1,4],1), volt([1,4]), 1e-14);
   tst = [ ...
   1.503131926779798; 1.412534629974291; 1.529078332819747;
   1.354399248512161; 1.546241676995996];
   unit_test_cmp('a2c2 Vstim #3', vh.volt(1:5:25,1), tst, 1e-14);

   img.fwd_model.stimulation(2) = stimulv;
   img.fwd_model.stimulation(1).stim_pattern([1,4]) = 0;
   vh = fwd_solve_1st_order(img);
   unit_test_cmp('a2c2 Vstim #4', vh.volt(1:5:25,2), tst, 1e-14);

   imgn = rmfield(img,'elem_data'); imgn.node_data = vh.volt;
   imgn.calc_colours.clim = 1; subplot(221); show_fem(imgn,1);

   img = mk_image( mk_common_model('a2C2',16),1);
   img.fwd_model.stimulation = stimulv;
   img.fwd_solve.get_all_nodes = 1;
   vh = fwd_solve_1st_order(img);
   unit_test_cmp('a2C2 Vstim #1', vh.meas, diff(volt(1:2)), 1e-14);
   unit_test_cmp('a2C2 Vstim #2', vh.volt(num_nodes(img)+[1,4]), volt([1,4]), 1e-14);
   tst = [ ...
   1.499999999999998; 1.302478674263331; 1.609665333411830; ...
   1.215039511028270; 1.691145536046686];
   unit_test_cmp('a2C2 Vstim #3', vh.volt(1:5:25), tst, 1e-13);

   imgn = rmfield(img,'elem_data'); imgn.node_data = vh.volt(1:num_nodes(img));
   imgn.calc_colours.clim = 1; subplot(222); show_fem(imgn,1);

   stim = zeros(16,1); volt=stim; stim([1,4]) = NaN; stim(8)=1; volt([1,4]) = [1,1];
   img.fwd_model.stimulation(2).stim_pattern = stim;
   img.fwd_model.stimulation(2).volt_pattern = volt;
   img.fwd_model.stimulation(2).meas_pattern = [1,-1,zeros(1,14)];
   vh = fwd_solve_1st_order(img);
   unit_test_cmp('a2C2 Vstim #4', vh.volt(num_nodes(img)+[1,4],1), [1;2], 1e-14);
   unit_test_cmp('a2C2 Vstim #5', vh.volt(1:5:25,1), tst, 1e-13);
   unit_test_cmp('a2C2 Vstim #6', vh.volt(num_nodes(img)+[1,4],2), [1;1], 1e-14);
   tst = [ 1.029942389400905; 1.024198991581187; ...
           1.048244746016660; 1.006551737030278; 1.057453501332724];
   unit_test_cmp('a2C2 Vstim #7', vh.volt(1:5:25,2), tst, 1e-13); % needs weaker tolerance

   imgn = rmfield(img,'elem_data'); imgn.node_data = vh.volt(1:num_nodes(img),2);
   subplot(223); show_fem(imgn,1);

   stim = zeros(16,1); volt=stim; stim([3,6]) = NaN;  volt([3,6]) = [1,2];
   img.fwd_model.stimulation(3).stim_pattern = stim;
   img.fwd_model.stimulation(3).volt_pattern = volt;
   img.fwd_model.stimulation(3).meas_pattern = [1,-1,zeros(1,14)];
   vh = fwd_solve_1st_order(img);

   imgn = rmfield(img,'elem_data'); imgn.node_data = vh.volt(1:num_nodes(img),3);
   imgn.calc_colours.clim = 1; subplot(224); show_fem(imgn,1);

   unit_test_cmp('a2C2 Vstim #7', vh.volt(num_nodes(img)+[3,6],3), [1;2], 1e-14);



function do_unit_test
   unit_test_voltage_stims;


   img = mk_image( mk_common_model('b2c2',16),1);
   vh = fwd_solve_1st_order(img);
   tst = [ ...
    0.959567140078593; 0.422175237237900; 0.252450963869202; ...
    0.180376116490602; 0.143799778367518];
   unit_test_cmp('b2c2 TEST', vh.meas(1:5), tst, 1e-12);

   img.fwd_model = rmfield(img.fwd_model,'gnd_node');
   vh = fwd_solve_1st_order(img);
   unit_test_cmp('b2c2 gnd_node', vh.meas(1:5), tst, 1e-12);

   img.fwd_solve.get_elec_curr = 1;
   vh = fwd_solve_1st_order(img);
   pp = fwd_model_parameters( img.fwd_model); EC = pp.N2E*pp.QQ;
   unit_test_cmp('b2b2 (CEM) elec_curr', vh.elec_curr, EC, 1e-11);

   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve_1st_order(img);
    plot(vh.volt);

   img = mk_image( mk_common_model('b2C2',16),1);
   vh = fwd_solve_1st_order(img);
   tst = [ 0.385629619754662; 0.235061644846908; 0.172837756982388
           0.142197580506776; 0.126808900182258; 0.120605655110661];
   unit_test_cmp('b2C2 (CEM) TEST', vh.meas(15:20), tst, 1e-12);

   img.fwd_solve.get_elec_curr = 1;
   vh = fwd_solve_1st_order(img);
   pp = fwd_model_parameters( img.fwd_model); EC = pp.N2E*pp.QQ;
   unit_test_cmp('b2C2 (CEM) elec_curr', vh.elec_curr, EC, 1e-11);

   % bad stim patterns (flow through ground node)
   img.fwd_model.stimulation(1).stim_pattern(2) = 0;
   vh = fwd_solve_1st_order(img);
   last_msg = eidors_msg('last_msg');
   last_msg = last_msg(end-68:end);
   unit_test_cmp('gnd_node warning', last_msg, ...
    ' of current is flowing through ground node. Check stimulation pattern');


   %2D resistor
   current = 4; measure=1;
   [R,img] = test_2d_resistor(current,measure);
   img.fwd_solve.get_all_nodes = 1;
   vs = fwd_solve_1st_order( img);
   va= measure*current*sum(R); % analytic
   unit_test_cmp('2D resistor test', va, vs.meas, 1e-12);

   unit_test_cmp('2D R z_contact', ...
                 [diff(vs.volt([13,1])), diff(vs.volt([14,12]))], ...
                 R(2)/2*current*[1,-1], 1e-12);
   unit_test_cmp('2D R voltages', vs.volt(1:3:10)-vs.volt(1), ...
                 R(1)*current*linspace(0,1,4)', 1e-12);

   [R,img] = test_2d_resistor_faces(current,measure);
   vs = fwd_solve_1st_order( img);
   unit_test_cmp('2D resistor faces', va, vs.meas, 1e-12);

   %3D resistor
   [R,img] = test_3d_resistor(current,measure);
   img.fwd_solve.get_all_nodes = 1;
   vs = fwd_solve_1st_order( img);
   va= current*sum(R);
   unit_test_cmp('3D resistor test', va, vs.meas, 1e-10);
   unit_test_cmp('3D R voltages', vs.volt(1:12:72)-vs.volt(1), ...
                 R(1)*current*linspace(0,1,6)', 1e-10);
   unit_test_cmp('3D R z_contact', ...
                 [diff(vs.volt([73,1])), diff(vs.volt([74,72]))], ...
                 R(2)/2*current*[1,-1], 1e-10);

   [R,img] = test_3d_resistor_faces(current,measure);
   vs = fwd_solve_1st_order( img);
   unit_test_cmp('3D resistor faces', va, vs.meas, 1e-10);


function [R,img] = test_2d_resistor(current,measure)
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 3; len = 12; 

   fmdl=mk_grid_model([],linspace(0,wid,3), linspace(0,len,4));
   fmdl.electrode(1).nodes = find(fmdl.nodes(:,2) ==   0);
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,2) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R = len / wid / conduc;
   Contact_R = z_contact/wid;
   R = [Block_R, 2*Contact_R];

% define electrode using face rather than nodes
function [R,img] = test_2d_resistor_faces(current,measure)
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 3; len = 12; 

   fmdl=mk_grid_model([],linspace(0,wid,3), linspace(0,len,4));
   bdy = fmdl.boundary;
   bdy( any(reshape(fmdl.nodes(bdy,2),size(bdy))>0,2),:)=[];
   fmdl.electrode(1).nodes = [];
   fmdl.electrode(1).faces = bdy;
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,2) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R = len / wid / conduc;
   Contact_R = z_contact/wid;
   R = [Block_R, 2*Contact_R];

function [R,img] = test_3d_resistor(current,measure);;
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 2; len = 5; hig=3; 

   fmdl=mk_grid_model([],0:wid, 0:hig, 0:len);
   fmdl.electrode(1).nodes = find(fmdl.nodes(:,3) ==   0);
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,3) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R =  len / wid / hig / conduc;
   Contact_R = z_contact/(wid*hig);
   R = [Block_R, 2*Contact_R];

% define electrode using face rather than nodes
function [R,img] = test_3d_resistor_faces(current,measure);;
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= .1; wid = 2; len = 5; hig=3; 

   fmdl=mk_grid_model([],0:wid, 0:hig, 0:len);
%  fmdl.electrode(1).nodes = find(fmdl.nodes(:,3) ==   0);
   bdy = fmdl.boundary;
   bdy( any(reshape(fmdl.nodes(bdy,3),size(bdy))>0,2),:)=[];
   fmdl.electrode(1).nodes = [];
   fmdl.electrode(1).faces = bdy;
   fmdl.electrode(2).nodes = find(fmdl.nodes(:,3) == len);
   [fmdl.electrode(:).z_contact] = deal(z_contact);
   fmdl.stimulation = stim_meas_list([1,2,1,2],2,current,measure);
   img= mk_image(fmdl,conduc);

   Block_R =  len / wid / hig / conduc;
   Contact_R = z_contact/(wid*hig);
   R = [Block_R, 2*Contact_R];

function [R,img] = test_3d_resistor_old(current,measure);
   ll=5*1; % length
   ww=1*2; % width
   hh=1*3; % height
   conduc= .13;  % conductivity in Ohm-meters
   z_contact= 1e-1;
   scale = .46;
   nn=0;
   for z=0:ll; for x=0:ww; for y=0:hh
      nn=nn+1;
      mdl.nodes(nn,:) = [x,y,z];
   end; end; end
   mdl= eidors_obj('fwd_model','3D rectangle');
   mdl= mk_grid_model([],0:ww,0:hh,0:ll);
   mdl.nodes= mdl.nodes*scale;
   mdl= rmfield(mdl,'coarse2fine');

   mdl.boundary= find_boundary(mdl.elems);
   mdl.gnd_node = 1;
   elec_nodes= [1:(ww+1)*(hh+1)];
   elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
   elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
   stim.stim_pattern= [-1;1]*current;
   stim.meas_pattern= [-1,1]*measure;
   mdl.stimulation= stim;
   mdl.electrode= elec;
   mdl = mdl_normalize(mdl,0);

   mdl.solve = @fwd_solve_1st_order;
   mdl.system_mat = @system_mat_1st_order;
   img= eidors_obj('image','3D rectangle', ...
         'elem_data', ones(size(mdl.elems,1),1) * conduc, ...
         'fwd_model', mdl); 

   % analytical solution
   Block_R =  ll / ww / hh / scale/ conduc;
   Contact_R = z_contact/(ww*hh)/scale^2;
   % Contact R reflects z_contact / (width/scale)^2. Here we need to use
   %  the scale, since this is not reflected in the size of the
   %  FEM as created by the grid. This is different to the test_2d_resistor,
   %  where the FEM is created scaled, so that the ww
   %  don't need to be scaled by the scale parameter.
   R =[ Block_R , 2*Contact_R];
