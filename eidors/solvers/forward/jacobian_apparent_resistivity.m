function J =jacobian_apparent_resistivity(fwd_model, img)
% jacobian_apparent_resistivity: Jacobian output as apparent resistivity
%  This function is a wrapper to the calc_jacobian; however, the output
%  is converted into apparent resistivity units, rather than in
%  voltage units

% $Id$

% correct input paralemeters if function was called with only img
if ischar(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_APPARENT_RESISTIVITY with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

if mdl_normalize(fwd_model)
   error('Cannot calculate apparent resistivity for normalized difference data');
end


solver = @eidors_default;
try
   solver = fwd_model.jacobian_apparent_resistivity.jacobian;
end

img.fwd_model.jacobian = solver;

J = calc_jacobian(img);

fctr = get_factor(img);

J = fctr * J;

function fctr = get_factor(img)
fctr = NaN;
try
   fctr = img.fwd_model.apparent_resistivity_factor;
end

if ischar(fctr), fctr = str2func(fcstr); end;

if isa(fctr, 'function_handle')
   fctr = feval(fctr, img);
end

if isnan(fctr)
   vh = fwd_solve(mk_image(img.fwd_model,1));
   n = length(vh.meas);
   fctr = spdiags(1./vh.meas,0,n,n);
end


function do_unit_test
   disp('TODO: NEED TO WRITE TESTS');
