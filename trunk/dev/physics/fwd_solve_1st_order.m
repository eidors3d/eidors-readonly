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

try 
    measurement = fwd_model.measured_quantity;
catch
    measurement = 'voltage';
end

if ~ismember(measurement, supported_measurement)
    error('EIDORS:MeasurementNotSupported', '%s does not support %s',...
    'FWD_SOLVE_1ST_ORDER',measurement);
end

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

vv = extract_volts( fwd_model.stimulation, pp, v_els);

if ~strcmp(measurement, 'voltage')
    [vv,fctr] = convert_measurement(vv,img,measurement);
    if ~isempty(fctr)
        data.apparent_resistivity_factor= fctr;
    end
end
    

% create a data structure to return
data.meas= vv;

data.measured_quantity = measurement;
data.time= NaN; % unknown
data.name= 'solved by fwd_solve_1st_order';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = v;                % all, including CEM nodes
end; end

function [vv,fctr] = convert_measurement(vv,img,measurement)
switch measurement
    case 'abs_voltage'
        vv = abs(vv);
        fctr = [];
    case 'log_voltage'
        vv = log(vv);
        fctr= [];
    case 'log10_voltage'
        vv = log10(vv);
        fctr= [];
    case 'apparent_resistivity'
        fctr = apparent_resistivity_factor(img.fwd_model);
        vv = fctr*vv;
    case 'log_apparent_resistivity'
        fctr = apparent_resistivity_factor(img.fwd_model);
        vv = log(fctr*vv);
    case 'log10_apparent_resistivity'
        fctr = apparent_resistivity_factor(img.fwd_model);
        vv = log10(fctr*vv);
end

function vv = extract_volts( stim, pp, v_els)
try
    if 0
       % measured voltages from v
       vv = zeros( pp.n_meas, 1 );
       idx=0;
       for i=1:pp.n_stim
          meas_pat= stim(i).meas_pattern;
          n_meas  = size(meas_pat,1);
          vv( idx+(1:n_meas) ) = meas_pat*v_els(:,i);
          idx= idx+ n_meas;
       end
    else
       v_out = arrayfun(@(s,i) ...
               s.meas_pattern * v_els(:,i), ...
            stim, 1:pp.n_stim, 'UniformOutput', false);
       vv = vertcat( v_out{:} );
    end
catch err
   if strcmp(err.identifier, 'MATLAB:innerdim');
       error(['measurement pattern not compatible with number' ...
               'of electrodes for stimulation patter %d'],i);
   else
       rethrow(err);
   end
end

function list = supported_measurement
   list = {'voltage'
           'abs_voltage'
           'log_voltage'
           'log10_voltage'
           'apparent_resistivity'
           'log_apparent_resistivity'
           'log10_apparent_resistivity'
           };

function do_unit_test
   img = mk_image( mk_common_model('b2c2',16),1);
   vho = fwd_solve_1st_order(img);
   tst= [ 0.959567140078593; 0.422175237237900; 0.252450963869202];
   unit_test_cmp('values',vho.meas(1:3),tst,1e-13);

   img.elem_data = 0.1 + rand(size(img.elem_data));
   vh = fwd_solve_1st_order(img);
   img.fwd_model.measured_quantity = 'log_voltage';
   vl = fwd_solve_1st_order(img);
   subplot(221)
   plot(log(vh.meas),vl.meas,'.'); 
   unit_test_cmp('log_voltage',log(vh.meas),vl.meas);
   
   img.fwd_model.measured_quantity = 'apparent_resistivity';
   vr = fwd_solve_1st_order(img);
   subplot(222)
   plot(vh.meas./vho.meas,vr.meas,'.'); 
   unit_test_cmp('apparent_resistivity',vh.meas./vho.meas,vr.meas,1e-5);
   
   img.fwd_model.measured_quantity = 'log_apparent_resistivity';
   vlr = fwd_solve_1st_order(img);
   subplot(223)
   plot(log(vr.meas),vlr.meas,'.'); 
   unit_test_cmp('log_apparent_resistivity',log(vr.meas),vlr.meas);
   
   img.fwd_model.measured_quantity = 'log10_apparent_resistivity';
   vlr = fwd_solve_1st_order(img);
   subplot(224)
   plot(log10(vr.meas),vlr.meas,'.'); 
   unit_test_cmp('log10_apparent_resistivity',log10(vr.meas),vlr.meas);
  
  
