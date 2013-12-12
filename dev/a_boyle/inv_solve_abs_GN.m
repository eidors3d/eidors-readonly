function img= inv_solve_abs_GN( inv_model, data1);
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
% (C) 2013 Alistair Boyle
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); img = do_unit_test; return; end

% fixed working data... otherwise we wouldn't be calling this function!
inv_model.parameters.elem_working = 'log_conductivity';
img = inv_solve_abs_core(inv_model, data1);

function pass = do_unit_test()
   pass = inv_solve_abs_core('UNIT_TEST');
