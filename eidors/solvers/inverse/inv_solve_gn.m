function img= inv_solve_gn( inv_model, data1, data2);
%function img= inv_solve_gn( inv_model, data1);
% INV_SOLVE_GN
% This function calls INV_SOLVE_CORE to find a Gauss-Newton
% iterative solution.
%
% img = inv_solve_gn( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => measurements
%
% Example:
%  imdl=mk_common_model('a2c');
%  imdl.reconst_type='absolute'; % ***
%  imdl.solve=@inv_solve_gn; % ***
%  fimg=mk_image(imdl,1);
%  fimg.elem_data(5:10)=1.1;
%  vi=fwd_solve(fimg);
%  img=inv_solve(imdl,vi); % ***
%  show_fem(img,1);
%
% See INV_SOLVE_CORE for arguments, options and parameters.
%
% (C) 2014 Alistair Boyle
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); img = do_unit_test; return; end


% inv_model.inv_solve_gn -> inv_solve_core
if isfield(inv_model, 'inv_solve_gn')
   if isfield(inv_model, 'inv_solve_core')
      error('inv_model.inv_solve_gn replaces inv_model.inv_solve_core, parameters will be lost');
   end
   inv_model.inv_solve_core = inv_model.inv_solve_gn;
   inv_model = rmfield(inv_model, 'inv_solve_gn');
end

if nargin > 2
  img = inv_solve_core(inv_model, data1, data2);
else
  img = inv_solve_core(inv_model, data1);
end

if isfield(img, 'inv_solve_core')
  img.inv_solve_gn = img.inv_solve_core;
  img=rmfield(img, 'inv_solve_core');
end

function pass = do_unit_test()
   pass=1;
   imdl=mk_common_model('a2c');
   imdl.reconst_type='absolute'; % ***
   imdl.solve=@inv_solve_gn; % ***
   fimg=mk_image(imdl,1);
   fimg.elem_data(5:10)=1.1;
   vi=fwd_solve(fimg);
   img=inv_solve(imdl,vi); % ***
   clf; subplot(121); show_fem(fimg,1); title('forward model');
        subplot(122); show_fem(img,1);  title('reconstruction');
   try unit_test_cmp('fwd vs. reconst', fimg.elem_data, img.elem_data, 0.08);
   catch me; disp(me.message); pass=0; end
%   pass = pass & inv_solve_core('UNIT_TEST', 'inv_solve_gn');
