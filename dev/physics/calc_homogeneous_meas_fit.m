%function bg = calc_homogeneous_meas_fit(fmdl, va)
%function bg = calc_homogeneous_meas_fit(fmdl, va, elem_type)
%
% Generate a homogeneous background estimate based on the measured values.
%   fmdl - a forward model
%   va   - measurements
%   elem_type - the output element type; conductivity, log_resistivity, etc.
%
% See also MAP_MEAS, convert_img_units.
%
% (C) 2014 Alistair Boyle
% Licenced under GPL version 2 or 3
% $Id$
function bg = calc_homogeneous_meas_fit(fmdl, va, elem_type)
  if nargin < 3
    elem_type = 'conductivity';
  end
  if ~isstruct(va) % 'va' should have been a struct
    q = fmdl.measured_quantity;
    va = map_meas(va, q, q);
  end
  eidors_msg('est. background %s', elem_type, 2);
  cache_obj = { fmdl, va };
  bg = eidors_obj('get-cache', cache_obj, 'calc_homogeneous_fit');
  if isempty(bg)
    eidors_msg('  ... cache miss', 2);
    bg = conductivity_fit(fmdl, va, elem_type);
    eidors_obj('set-cache', cache_obj, 'calc_homogeneous_fit', bg);
  else
    eidors_msg('  ... cache hit', 2);
  end

  if 1
    bg_tmp = convert_img_units(bg, elem_type, 'resistivity');
    eidors_msg('  estimated background resistivity: %0.1f Ohm.m', bg_tmp, 2);
  end

function bg = conductivity_fit(fmdl, va, elem_type)
  % have a look at what we've created
  % compare data to homgeneous (how good is the model?)
  % NOTE background conductivity is set by matching amplitude of
  % homogeneous data against the measurements to get a rough fit
  imgh = mk_image(fmdl, 1); % conductivity = 1 S/m
  vh = fwd_solve(imgh);

  % take the best fit of the data (linear fit!)
  %   pf = polyfit(data,vs.meas,1); % <-- Nolwenn's solution is pf(1)
  vh = map_meas(vh, va.measured_quantity);
  bg = vh.meas \ va.meas; % <-- Alistair's solution

  % convert to output units
  bg = convert_img_units(bg, 'conductivity', elem_type);
