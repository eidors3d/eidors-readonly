function mr = calc_model_reducton(fmdl)
% calc_model_reduction: calculate the fields for a reduced model
%    which should speed up forward calculations
% 
% mr = calc_model_reducton(fmdl)
%  where
%    mr.main_region = vector, and 
%    mr.regions = struct
%
% Model Reduction: use precomputed fields to reduce the size of
%    the forward solution. Nodes which are 1) not used in the output
%    (i.e. not electrodes) 2) all connected to the same conductivity via
%    the c2f mapping are applicable.

% see: Model Reduction for FEM Forward Solutions, Adler & Lionheart, EIT2016

% $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

switch fmdl(1).type
  case 'image';      fmdl = fmdl.fwd_model;
  case 'inv_model';  fmdl = fmdl.fwd_model;
  case 'fwd_model';  fmdl = fmdl;
  otherwise;
      error('can''t process model of type %s', fmdl.type );
end


function do_unit_test
