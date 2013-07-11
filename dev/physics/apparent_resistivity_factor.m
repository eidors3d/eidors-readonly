function fctr = apparent_resistivity_factor(fmdl)

switch fmdl.type
   case 'fwd_model' % ok
   case 'image' 
      fmdl = fmdl.fwd_model;
   otherwise
      error('EIDORS:WrongInput','Image or fwd_model struct required.');
end

fctr = eidors_cache(@calc_factor,{fmdl},'apparent_resistivity_factor');

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