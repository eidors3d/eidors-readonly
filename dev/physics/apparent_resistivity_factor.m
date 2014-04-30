function fctr = apparent_resistivity_factor(fmdl)

switch fmdl.type
   case 'fwd_model' % ok
   case 'image' 
      fmdl = fmdl.fwd_model;
   otherwise
      error('EIDORS:WrongInput','Image or fwd_model struct required.');
end

if ~isfield(fmdl,'apparent_resistivity_factor');
  fmdl.apparent_resistivity_factor = NaN;
end

if ~isfield(fmdl, 'normalize_measurements')
  fmdl = mdl_normalize(fmdl,0);
end

fmdlfields = fieldnames(fmdl);
keepfields = {'type','elems','nodes','solve','system_mat','gnd_node', ...
  'electrode','stimulation','normalize_measurements','apparent_resistivity_factor'};

keep = ismember(fmdlfields, keepfields);

%only cache on these fields
cache_obj = rmfield(fmdl, fmdlfields(~keep));

fctr = eidors_cache(@calc_factor,{cache_obj},'apparent_resistivity_factor');

function fctr = calc_factor(fmdl)

fctr = NaN;
try
   fctr = fmdl.apparent_resistivity_factor;
end

if ischar(fctr), fctr = str2func(fcstr); end;

if isa(fctr, 'function_handle')
   fctr = feval(fctr, img);
end

if isnan(fctr)
   fmdl.measured_quantity = 'voltage';
   vh = fwd_solve(mk_image(fmdl,1));
   n = length(vh.meas);
   fctr = spdiags(1./vh.meas,0,n,n);
end