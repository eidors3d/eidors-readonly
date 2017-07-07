function img= inv_solve_cg( inv_model, data1, data2);
%function img= inv_solve_cg( inv_model, data1);
% INV_SOLVE_CG
% This function calls INV_SOLVE_ABS_CORE to find a Conjugate-Gradient
% iterative solution.
%
% img = inv_solve_cg( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => measurements
%
% Example:
%  imdl=mk_common_model('a2c');
%  imdl.reconst_type='absolute'; % ***
%  imdl.solve=@inv_solve_cg; % ***
%  fimg=mk_image(imdl,1);
%  fimg.elem_data(5:10)=1.1;
%  vi=fwd_solve(fimg);
%  img=inv_solve(imdl,vi); % ***
%  show_fem(img,1);
%
% See INV_SOLVE_ABS_CORE for arguments, options and parameters.
%
% (C) 2014 Alistair Boyle
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

% fixed working data... otherwise we wouldn't be calling this function!
if ~isfield(inv_model, 'inv_solve_cg')
   inv_model.inv_solve_cg = struct;
end
if ~isfield(inv_model.inv_solve_cg, 'update_func')
   inv_model.inv_solve_cg.update_func = 'GN_update';
end
if ~isfield(inv_model.inv_solve_cg, 'beta_func')
   inv_model.inv_solve_cg.beta_func = @beta_reset_polak_ribiere;
end

% inv_model.inv_solve_cg -> inv_solve_core
if isfield(inv_model, 'inv_solve_cg')
   inv_model.inv_solve_core = inv_model.inv_solve_cg;
   inv_model = rmfield(inv_model, 'inv_solve_cg');
end

if nargin > 2
   img = inv_solve_core(inv_model, data1, data2);
else
   img = inv_solve_core(inv_model, data1);
end

if isfield(img, 'inv_solve_core')
  img.inv_solve_cg = img.inv_solve_core;
  img=rmfield(img, 'inv_solve_core');
end

% Note: we prefer to still be able to apply regularization, so use the
% GN_update by default, but you can try this if you are determined. It
% generally doesn't give good reconstructions for most EIT inverse problems.
function dx = CG_update(J, W, hp2RtR, dv, de, opt)
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

% Dai-Yuan
% cite: http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
%                dx_k^T    dx_k
% beta = - -------------------------------
%          s_{k-1}^T ( dx_k - dx_{k-1} )
function beta = beta_dai_yuan(dx_k, dx_km1, sx_km1)
   beta = -dx_k'*dx_k/(sx_km1'*(dx_k-dx_km1));

function do_unit_test()
   imdl=mk_common_model('a2c');
   imdl.reconst_type='absolute'; % ***
   imdl.solve=@inv_solve_cg; % ***
   fimg=mk_image(imdl,1);
   fimg.elem_data(5:10)=1.1;
   vi=fwd_solve(fimg);
   img=inv_solve(imdl,vi); % ***
   clf; subplot(121); show_fem(fimg,1); title('forward model');
        subplot(122); show_fem(img,1);  title('reconstruction');
   unit_test_cmp('fwd vs. reconst', fimg.elem_data, img.elem_data, 0.08);
%  inv_solve_core('UNIT_TEST', 'inv_solve_cg');
