function data =fwd_solve_apparent_resistivity(fwd_model, img)
% fwd_solve_apparent_resistivity: fwd_solve output as apparent resistivity
%  This function is a wrapper to the fwd_solve; however, the output
%  is converted into apparent resistivity units, rather than in
%  voltage units

% $Id$

% correct input paralemeters if function was called with only img

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling FWD_SOLVE_APPARENT_RESISTIVITY with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

solver = @eidors_default;
try
   solver = fwd_model.fwd_solve_apparent_resistivity.solve;
end

img.fwd_model.solve = solver;

data = fwd_solve(img);

fctr = get_factor(img);

data.meas = fctr * data.meas;
data.name = ['apparent resistivity ' data.name];
data.quantity = 'apparent resistivity';
data.apparent_resistivity_factor= fctr;

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
   imdl = mk_common_model('a2c2',8); 
   img = mk_image(imdl);
   vrh = fwd_solve( img );
   img.elem_data(1) = 1.1;
   vri = fwd_solve( img );

   img = mk_image(imdl);
   img.fwd_model.solve = @fwd_solve_apparent_resistivity;
   vah = fwd_solve( img );
   img.elem_data(1) = 1.1;
   vai = fwd_solve( img );

   unit_test_cmp('homog is ones', vah.meas, ones(size(vah.meas)), 1e-10);
   unit_test_cmp('ratio is same', vai.meas./vah.meas, vri.meas./vrh.meas, 1e-10);
