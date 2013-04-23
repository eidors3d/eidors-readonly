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
if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling FWD_SOLVE_1ST_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

img = physics_data_mapper(img);
if ~ismember(img.current_physics, supported_physics)
    error('EIDORS:PhysicsNotSupported', '%s does not support %s', ...
    'FWD_SOLVE_1ST_ORDER',img.current_physics);
end
orig_physics = img.current_physics;
% all calcs use conductivity
img = convert_units(img, 'conductivity');

pp= fwd_model_parameters( fwd_model );
s_mat= calc_system_mat( img );

idx= 1:size(s_mat.E,1);
idx( fwd_model.gnd_node ) = [];

v= zeros(pp.n_node,pp.n_stim);

v(idx,:)= left_divide( s_mat.E(idx,idx), pp.QQ(idx,:));

% calc voltage on electrodes
v_els= pp.N2E * v;

try
    % measured voltages from v
    vv = zeros( pp.n_meas, 1 );
    idx=0;
    for i=1:pp.n_stim
       meas_pat= fwd_model.stimulation(i).meas_pattern;
       n_meas  = size(meas_pat,1);
       vv( idx+(1:n_meas) ) = meas_pat*v_els(:,i);
       idx= idx+ n_meas;
    end
catch err
   if strcmp(err.identifier, 'MATLAB:innerdim');
       error(['measurement pattern not compatible with number' ...
               'of electrodes for stimulation patter %d'],i);
   else
       rethrow(err);
   end
end


% create a data structure to return
data.meas= vv;
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_1st_order';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = v;                % all, including CEM nodes
end; end

function do_unit_test
   img = mk_image( mk_common_model('b2c2',16),1);
   vh = fwd_solve_1st_order(img);
