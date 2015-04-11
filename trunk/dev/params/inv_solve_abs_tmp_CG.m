function img= inv_solve_abs_cg( inv_model, data1);
%function img= inv_solve_abs_cg( inv_model, data1);
% INV_SOLVE_ABS_CG
% This function calls INV_SOLVE_ABS_CORE to find a Conjugate-Gradient
% iterative solution.
%
% img = inv_solve_abs_cg( inv_model, data1 )
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

% fixed working data... otherwise we wouldn't be calling this function!
if ~isfield(inv_model.inv_solve, 'update_func')
  inv_model.inv_solve.update_func = @CG_update;
end
if ~isfield(inv_model.inv_solve, 'beta_func')
  inv_model.inv_solve.beta_func = @beta_reset_polak_ribiere;
end

img = inv_solve_abs_core(inv_model, data1);

function dx = CG_update(J, W, hp2, RtR, dv, de)
%  hp2RtR = hp2*RtR;
%  % the actual update
%  dx = (J'*W*J + hp2RtR)\(J'*dv + hp2RtR*de);
  dx = J'*dv;

% Fletcher-Reeves
% cite: http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
%            dx_k^T  dx_k
% beta = ---------------------
%         dx_{k-1}^T dx_{k-1}
function beta = beta_fletcher_reeves(dx_k, dx_km1, sx_km1)
  beta = (dx_k' * dx_k)/(dx_km1' * dx_km1);
  if isinf(beta) || isnan(beta)
    beta = 0;
  end

% Polak-Ribiere
% cite: http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
%         dx_k^T ( dx_k - dx_{k-1} )
% beta = ----------------------------
%             dx_{k-1}^T dx_{k-1}
function beta = beta_polak_ribiere(dx_k, dx_km1, sx_km1)
  ddx = dx_k - dx_km1;
  beta = (dx_k' * ddx)/(dx_km1' * dx_km1);
  if isinf(beta) || isnan(beta)
    beta = 0;
  end

% Hestenes-Stiefel
% cite: http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
%          dx_k^T    ( dx_k - dx_{k-1} )
% beta = - -------------------------------
%          s_{k-1}^T ( dx_k - dx_{k-1} )
function beta = beta_hestenes_stiefel(dx_k, dx_km1, sx_km1)
  ddx = dx_k - dx_km1;
  beta = -(dx_k' * ddx)/(sx_km1' * ddx);
  if isinf(beta) || isnan(beta)
    beta = 0;
  end

% Polak-Ribiere with reset
% cite: http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
% beta = max(0, beta_pr)
function beta = beta_reset_polak_ribiere(dx_k, dx_km1, sx_km1)
  beta = beta_polak_ribiere(dx_k, dx_km1, sx_km1);
  if beta < 0
    beta = 0;
  end

function pass = do_unit_test()
   pass = inv_solve_abs_core('UNIT_TEST', 'inv_solve_abs_cg');
