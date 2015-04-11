function img= inv_solve_abs_GN( inv_model, data1);
%function img= inv_solve_abs_GN( inv_model, data1);
% INV_SOLVE_ABS_GN
% This function calls INV_SOLVE_ABS_CORE to find a Gauss-Newton
% iterative solution.
%
% img = inv_solve_abs_GN( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => measurements
%
% See INV_SOLVE_ABS_CORE for arguments, options and parameters.
%
% (C) 2014 Alistair Boyle
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); img = do_unit_test; return; end

img = inv_solve_abs_core(inv_model, data1);

function pass = do_unit_test()
   pass = inv_solve_abs_core('UNIT_TEST', 'inv_solve_abs_GN');
