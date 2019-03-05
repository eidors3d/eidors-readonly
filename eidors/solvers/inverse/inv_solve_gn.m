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
% Example#2:
%imdl = mk_common_model('b2c2',16);
%imdl.jacobian_bkgnd.value = 2;
%img=mk_image(imdl);    vh=fwd_solve(img);
%img.elem_data([5])=5;   vi=fwd_solve(img);
%imgr = inv_solve(imdl,vh,vi);  % regular solve
%
%imdl.solve = @inv_solve_gn;
%imdl.inv_solve_gn.max_iterations = 5;
%imdl.inv_solve_gn.min_value = 2;
%imgr = inv_solve(imdl,vh,vi); % new solve
%
% See INV_SOLVE_CORE for arguments, options and parameters.
%
% (C) 2014 Alistair Boyle
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end


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

function do_unit_test()
   imdl=mk_common_model('a2c');
   imdl.reconst_type='absolute'; % ***
   imdl.solve=@inv_solve_gn; % ***
   fimg=mk_image(imdl,1);
   fimg.elem_data(5:10)=1.1;
   vi=fwd_solve(fimg);
   img=inv_solve(imdl,vi); % ***
   clf; subplot(121); show_fem(fimg,1); title('forward model');
        subplot(122); show_fem(img,1);  title('reconstruction');
   unit_test_cmp('fwd vs. reconst', fimg.elem_data, img.elem_data, 0.08);

   do_unit_test_abs_diff();

function do_unit_test_abs_diff()
   imdl = mk_geophysics_model('h2c',32);
   imdl.solve = 'inv_solve_gn';
   imdl.RtR_prior = 'prior_tikhonov';
   imdl.inv_solve_gn.verbose = 0;
   imdl.inv_solve_gn.elem_working = 'log_conductivity';
%   imdl.inv_solve_gn.meas_working = 'apparent_resistivity';
   % build fwd simulations
   ctrs = interp_mesh(imdl.fwd_model);
   x = ctrs(:,1);
   y = ctrs(:,2);
   r1 = sqrt((x-20).^2  + (y+25).^2);
   r2 = sqrt((x-125).^2 + (y+45).^2);
   r3 = sqrt((x-50).^2  + (y+35).^2);
   imga = mk_image(imdl,1);
   imga.elem_data(r1<15)= 0.5;
   imga.elem_data(r2<20)= 2;
   va = fwd_solve(imga);
   imgb = imga;
   imgb.elem_data(r3<15)= 1.3;
   vb = fwd_solve(imgb);
   clf; subplot(121); show_fem(imga,1); subplot(122); show_fem(imgb,1); drawnow;
   % reconstruct meas va as absolute
   disp('*** INVERSE CRIME ALERT **** INVERSE CRIME ALERT ***');
   disp(' 1. This reconstruction has no measurement noise.');
   disp(' 2. This reconstruction is being performed on the same mesh the');
   disp('    measurements were simulated from: no discretization errors have');
   disp('    been introduced.');
   imdl.reconst_type = 'absolute';
   imgaa = inv_solve(imdl,va);
   % set background as abs sol'n, then diff solve
   imdl.reconst_type = 'difference';
   imdl.inv_solve_gn.prior_data = imgaa.elem_data;
   imgab = inv_solve(imdl,va,vb);
   clf; subplot(221); show_fem(imga,1); title('A'); subplot(222); show_fem(imgb,1); title('B');
        subplot(223); show_fem(imgaa,1); title('rec A'); subplot(224); show_fem(imgab,1); title('rec B-A');
   unit_test_cmp('abs ', imga.elem_data, imgaa.elem_data, max(imga.elem_data)/2);
   imgba = imga; imgba.elem_data = imgb.elem_data - imga.elem_data;
   imgba.elem_data = imgba.elem_data/max(imgba.elem_data);
   imgab.elem_data = imgab.elem_data/max(imgab.elem_data);
   unit_test_cmp('diff', norm(imgba.elem_data - imgab.elem_data), 0, 20);
