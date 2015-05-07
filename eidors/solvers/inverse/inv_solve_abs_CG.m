function img= inv_solve_abs_CG( inv_model, data1);
%function img= inv_solve_abs_CG( inv_model, data1);
% INV_SOLVE_ABS_CG
% This function calls INV_SOLVE_ABS_CORE to find a Conjugate-Gradient
% iterative solution.
%
% img = inv_solve_abs_CG( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => measurements
%
% Example:
%  imdl=mk_common_model('a2c');
%  imdl.reconst_type='absolute'; % ***
%  imdl.solve=@inv_solve_abs_CG; % ***
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
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); img = do_unit_test; return; end

% merge legacy options locations
inv_model = deprecate_imdl_opt(inv_model, 'parameters');
inv_model = deprecate_imdl_opt(inv_model, 'inv_solve');
inv_model = deprecate_imdl_opt(inv_model, 'inv_solve_abs_core');

% fixed working data... otherwise we wouldn't be calling this function!
if ~isfield(inv_model, 'inv_solve_abs_CG')
   inv_model.inv_solve_abs_CG = struct;
end
if ~isfield(inv_model.inv_solve_abs_CG, 'update_func')
   inv_model.inv_solve_abs_CG.update_func = @CG_update;
end
if ~isfield(inv_model.inv_solve_abs_CG, 'beta_func')
   inv_model.inv_solve_abs_CG.beta_func = @beta_reset_polak_ribiere;
end

% inv_model.inv_solve_abs_CG -> inv_solve_abs_core
if isfield(inv_model, 'inv_solve_abs_CG')
   inv_model.inv_solve_abs_core = inv_model.inv_solve_abs_CG;
   inv_model = rmfield(inv_model, 'inv_solve_abs_CG');
end

img = inv_solve_abs_core(inv_model, data1);

if isfield(img, 'inv_solve_abs_core')
  img.inv_solve_abs_CG = img.inv_solve_abs_core;
  img=rmfield(img, 'inv_solve_abs_core');
end

function dx = CG_update(J, W, hp2RtR, dv, de)
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

% Dai-Yuan
% cite: http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
%                dx_k^T    dx_k
% beta = - -------------------------------
%          s_{k-1}^T ( dx_k - dx_{k-1} )
function beta = beta_dai_yuan(dx_k, dx_km1, sx_km1)
   beta = -dx_k'*dx_k/(sx_km1'*(dx_k-dx_km1));

function imdl = deprecate_imdl_opt(imdl,opt)
   if ~isfield(imdl, opt)
      return;
   end
   if ~isstruct(imdl.(opt))
      error(['unexpected inv_model.' opt ' where ' opt ' is not a struct... i do not know what to do']);
   end

   % warn on anything but inv_model.inv_solve.calc_solution_error
   Af = fieldnames(imdl.(opt));
   if ~strcmp(opt, 'inv_solve') || (length(Af(:)) ~= 1) || ~strcmp(Af(:),'calc_solution_error')
      disp(imdl)
      disp(imdl.(opt))
      warning('EIDORS:deprecatedParameters',['INV_SOLVE inv_model.' opt '.* are deprecated in favor of inv_model.inv_solve_abs_CG.* as of 30-Apr-2014.']);
   end

   if ~isfield(imdl, 'inv_solve_abs_CG')
      imdl.inv_solve_abs_CG = imdl.(opt);
   else % we merge
      % merge struct trick from:
      %  http://stackoverflow.com/questions/38645
      for i = fieldnames(imdl.(opt))'
         imdl.inv_solve_abs_CG.(i{1})=imdl.(opt).(i{1});
      end
   end
   imdl = rmfield(imdl, opt);

function pass = do_unit_test()
   pass=1;
   imdl=mk_common_model('a2c');
   imdl.reconst_type='absolute'; % ***
   imdl.solve=@inv_solve_abs_CG; % ***
   fimg=mk_image(imdl,1);
   fimg.elem_data(5:10)=1.1;
   vi=fwd_solve(fimg);
   img=inv_solve(imdl,vi); % ***
   clf; subplot(121); show_fem(fimg,1); title('forward model');
        subplot(122); show_fem(img,1);  title('reconstruction');
   try unit_test_cmp('fwd vs. reconst', fimg.elem_data, img.elem_data, 0.08);
   catch me; disp(me.message); pass=0; end
%   pass = inv_solve_abs_core('UNIT_TEST', 'inv_solve_abs_CG');
