function img= inv_solve_abs_GN_logc( inv_model, data1);
% INV_SOLVE_ABS_GN_LOGC
% This function calls INV_SOLVE_ABS_CORE to find a Gauss-Newton
% iterative solution using log conductivity.
%
% img = inv_solve_abs_GN_logc( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => EIT measurements
%
% See INV_SOLVE_ABS_CORE for arguments, options and parameters.
%
% (C) 2013 Alistair Boyle, Nolwenn Lespare, Andy Adler
% License: GPL version 2 or version 3
% $Id$

% we assume resistivity output...
% probably should be conductivity to be consistent with the rest of EIDORS
%if ~isfield(inv_model.parameters, 'elem_output')
%   inv_model.parameters.elem_output = 'resistivity';
%end

% fixed working data... otherwise we wouldn't be calling this function!
inv_model.parameters.elem_working = 'log_conductivity';
img = inv_solve_abs_core(inv_model, data1);
