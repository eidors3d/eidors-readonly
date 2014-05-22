function test_params(N)
% Test the new parametrization interface

if nargin==0
   N = 1:12;
end
for i = N
   run_test(i);
end

function run_test(N)
eidors_msg('log_level',0);
eidors_cache off
imdl = prepare_model;
fmdl = imdl.fwd_model;
switch N
   case 1
      imgh = mk_image(fmdl,1,'Old code default');
      imgi = imgh;
      imgi.elem_data(32) = 2;
   case 2
      imgh = mk_image(fmdl,1,'conductivity','New code default');
      imgi = imgh;
      imgi.conductivity.elem_data(32) = 2;
   case 3
      imgh = mk_image(fmdl,1,'conductivity','Conductivity');
      imgi = imgh;
      imgi.conductivity.elem_data(32) = 2;
      imdl.jacobian_bkgnd = imgh;
   case 4
      imgh = mk_image(fmdl,1,'resistivity','Resistivity');
      imgi = imgh;
      imgi.resistivity.elem_data(32) = 1/2;
      imdl.jacobian_bkgnd = imgh;
   case 5
      fmdl.jacobian = @dummy_jacobian;
      fmdl.solve    = @dummy_fwd_solve;
      imgh = mk_image(fmdl,1,'resistivity','Custom functions');
      imgi = imgh;
      imgi.resistivity.elem_data(32) = 1/2;
      imdl.jacobian_bkgnd = imgh;
      imdl.fwd_model = fmdl;
   case 6
      imdl = mk_common_model('b3cr',[16,1]);
      fmdl = imdl.fwd_model;
      imgh = mk_image(fmdl,1,'resistivity','GREIT resistivity');
      imgi = imgh;
      imgi.resistivity.elem_data(3450) = 1/2;
      imdl = mk_GREIT_model(imgh,0.05,2);
   case 7
      imdl = mk_common_model('b3cr',[16,1]);
      fmdl = imdl.fwd_model;
      imgh = mk_image(fmdl,1,'resistivity','GREIT resistivity');
      imgi = imgh;
      imgi.resistivity.elem_data(3450) = 1/2;
      subplot(4,4,N)
      show_slices(imgi,[inf inf -0.1]);
      title('GREIT resistivity target');
      return
   case 8
      % graphics functions
      imdl = mk_common_model('b3cr',[16,1]);
      fmdl = imdl.fwd_model;
      imgh = mk_image(fmdl,1,'resistivity','Graphics functions');
      imgh.resistivity.elem_data(3400:3500) = 1/2;
      if isfield(imgh,'elem_data'); error('bad'); end
      subplot(4,4,N)
      show_fem(imgh)
      show_slices(imgh,[inf inf -0.1])
      show_3d_slices(imgh)
      title('Graphics functions');
      return
   case 9
      imdl.solve = @inv_solve_abs_GN;
      imdl.reconst_type = 'absolute';
      imdl.parameters.max_iterations = 3;
      imgh = mk_image(fmdl,1,'Absolute old');
      imgi = imgh;
      imgi.elem_data(32) = 2;
   case 10
      imdl.solve = @inv_solve_abs_GN;
      imdl.reconst_type = 'absolute';
      imdl.parameters.max_iterations = 3;
      imdl = rmfield(imdl,'jacobian_bkgnd');
      imdl.jacobian_bkgnd.conductivity.elem_data = 1;
      imgh = mk_image(fmdl,1,'Absolute new');
      imgi = imgh;
      imgi.elem_data(32) = 2;
   case 11
      imdl.solve = @inv_solve_abs_GN;
      imdl.reconst_type = 'absolute';
      imdl.parameters.max_iterations = 3;
      imdl = rmfield(imdl,'jacobian_bkgnd');
      imdl.jacobian_bkgnd.resistivity.elem_data = 1;
      imgh = mk_image(fmdl,1,'Absolute res');
      imgi = imgh;
      imgi.elem_data(32) = 2;
      % imgi is a conductivity image, just because it can
   case 12
      test_GREIT_resistivity;
      return
   otherwise
      error('No test %d',N);
end
fprintf([imgh.name '...\t']);
try
   %     S = calc_system_mat(imgh); %???
   J = calc_jacobian(imgh);
   vh = fwd_solve(imgh);
   vi = fwd_solve(imgi);
   if N < 9
      imgr = inv_solve(imdl,vh, vi);
      imgr.calc_colours.ref_level = 0;
      imgr.calc_colours.clim = 0.05;
   else
      imgr = inv_solve(imdl,vi);
   end
   fprintf('OK\n');
catch err
   if strcmp(err.identifier,'EIDORS:PhysicsNotSupported')
      fprintf('NOT SUPPORTED\n');
      return
   else
      rethrow(err)
   end
end
subplot(4,4,N)
show_slices(imgr);
eidors_colourbar(imgr);
title(imgh.name);
eidors_msg('log_level',2);
eidors_cache on


function [imgh imgi imdl] = test_GREIT_resistivity

x = 0.8;

imdl = mk_common_model('b3cr',[16,1]);
fmdl = imdl.fwd_model;
imgh = mk_image(fmdl,1,'resistivity','GREIT resistivity');
fm.elem_centre = 1;
fmdl = fix_model(fmdl,fm);
imgh.resistivity.elem_data(fmdl.elem_centre(:,2) > 0) = 1 + x;
imgh.resistivity.elem_data(fmdl.elem_centre(:,2) > 0) = 1 - x;
imgi = imgh;
select_fcn = inline('(x).^2+(y).^2+(z).^2<0.2^2','x','y','z');
memb_frac = elem_select( imgi.fwd_model, select_fcn);
imgi.resistivity.elem_data = imgi.resistivity.elem_data + memb_frac*1;
tmp = mk_image(fmdl,1,'resistivity','GREIT resistivity');
imdl = mk_GREIT_model(tmp,0.05,2);
% 
% imdl = mk_common_model('b3cr',[16,1]);
% imdl.jacobian_bkgnd = imgh;

vh = fwd_solve(imgh);
vi = fwd_solve(imgi);
imgr = inv_solve(imdl,vh, vi);
imgr.calc_colours.ref_level = 0;
subplot(4,4,12)
show_slices(imgr,[inf inf 0]);

function imdl = prepare_model
imdl = mk_common_model('b2C',16);



