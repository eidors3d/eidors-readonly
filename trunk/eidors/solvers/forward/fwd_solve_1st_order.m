function data =fwd_solve_1st_order(fwd_model, img)
% FWD_SOLVE_1ST_ORDER: data= fwd_solve_1st_order( img)
% Fwd solver for Andy Adler's EIT code
% Input:
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)

% (C) 1995-2002 Andy Adler. License: GPL version 2 or version 3
% Ref: Adler & Guardo (1996) IEEE T. Med Imaging
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
s_mat= calc_system_mat( img );

idx= 1:size(s_mat.E,1);
gnd_node = find_gnd_node( fwd_model );
idx( gnd_node ) = [];

v= zeros(pp.n_node,pp.n_stim);

v(idx,:)= left_divide( s_mat.E(idx,idx), pp.QQ(idx,:));

% calc voltage on electrodes

% This is horribly inefficient, override
% v_els= pp.N2E * v;
idx = find(any(pp.N2E));
v_els= pp.N2E(:,idx) * v(idx,:);


% create a data structure to return
data.meas= meas_from_v_els(v_els, fwd_model.stimulation);
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_1st_order';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = v;                % all, including CEM nodes
end; end


function gnd_node = find_gnd_node( fwd_model );
   if isfield(fwd_model,'gnd_node')
      gnd_node = fwd_model.gnd_node;
   else
      error('no gnd_node on fwd_model')
   end

function vv = meas_from_v_els( v_els, stim)
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

       copt.cache_obj = {stim};
       copt.fstr = 'v2meas';
       copt.log_level = 4;
       v2meas = eidors_cache(@get_v2meas, {n_elec,n_stim,stim}, copt);
       vv = v2meas' * v_els(:);
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
        v2meas((i-1)*n_elec + 1: i*n_elec,end+(1:n_meas)) = meas_pat';
    end
        

function do_unit_test
   img = mk_image( mk_common_model('b2c2',16),1);
   vh = fwd_solve_1st_order(img);
plot(vh.meas);
   unit_test_cmp('b2c2 TEST', vh.meas(1:5), ...
  [ 0.959567140078593; 0.422175237237900; 0.252450963869202; ...
    0.180376116490602; 0.143799778367518], 1e-8);
   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve_1st_order(img);
plot(vh.volt);

   test_2d_resistor;
   test_3d_resistor;

function test_2d_resistor
   current= 4;  % Amps
   conduc=  .4 + 2*pi*j*10; % conductivity in Ohm-meters
   z_contact= 1e-1;
   nn= 12;     % number of nodes
   ww=3;       % width = 4
   scale = .35;
   mdl=mk_grid_model([],3+scale*(1:ww), scale*(1:nn/ww));
   mdl= rmfield(mdl,'coarse2fine'); % don't calc this.

   mdl.gnd_node = 1;
   elec_nodes= [1:ww];
   elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
   elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
   stim.stim_pattern= [-1;1]*current;
   stim.meas_pattern= [-1,1];
   mdl.stimulation= stim;
   mdl.electrode= elec;
   n_el = size(mdl.elems,1);
   img= eidors_obj('image','2D rectangle', ...
         'elem_data', ones(n_el,1) * conduc );
   img.fwd_model = mdl;
   img.fwd_model.normalize_measurements = 0;
   img.fwd_model.solve = @fwd_solve_1st_order;
   img.fwd_model.system_mat = @system_mat_1st_order;

   vs = fwd_solve( img);

% Analytic
   nodes = img.fwd_model.nodes;
   wid_len= max(nodes) - min(nodes);
   conduc =  mean(img.elem_data);
   Block_R = wid_len(2) / wid_len(1) / conduc;
   % Contact R reflects z_contact / width. There is no need to scale
   %  by the scale, since this is already reflected in the size of the
   %  FEM as created by the grid. This is different to the test_3d_resistor,
   %  where the FEM is created first, and then scaled, so that the ww
   %  and hh need to be scaled by the scale parameter.
   z_contact = sum([img.fwd_model.electrode(:).z_contact]);
   current = max(img.fwd_model.stimulation(1).stim_pattern(:));
   Contact_R = z_contact/wid_len(1);
   R = Block_R + Contact_R;

   va= current*R;

   unit_test_cmp('2D resistor test', va, vs.meas, 1e-15);

function test_3d_resistor
   ll=5*1; % length
   ww=1*2; % width
   hh=1*3; % height
   conduc= .13;  % conductivity in Ohm-meters
   current= 4;  % Amps
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
   stim.meas_pattern= [-1,1];
   mdl.stimulation= stim;
   mdl.electrode= elec;
   mdl = mdl_normalize(mdl,0);

   mdl.solve = @fwd_solve_1st_order;
   mdl.system_mat = @system_mat_1st_order;
   img= eidors_obj('image','3D rectangle', ...
         'elem_data', ones(size(mdl.elems,1),1) * conduc, ...
         'fwd_model', mdl); 

   vs= fwd_solve(img);

   % analytical solution
   Block_R =  ll / ww / hh / scale/ conduc;
   Contact_R = z_contact/(ww*hh)/scale^2;
   % Contact R reflects z_contact / (width/scale)^2. Here we need to use
   %  the scale, since this is not reflected in the size of the
   %  FEM as created by the grid. This is different to the test_2d_resistor,
   %  where the FEM is created scaled, so that the ww
   %  don't need to be scaled by the scale parameter.
   R = Block_R + 2*Contact_R;
   va= current*R;
   unit_test_cmp('3D resistor test', va, vs.meas, 1e-10);
