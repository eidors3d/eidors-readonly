function test_physics(N)
% Test the new physics interface

if nargin==0
   N = 1:9;
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
      subplot(3,3,N)
      show_slices(imgi,[inf inf -0.1]);
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
   imgr = inv_solve(imdl,vh, vi);
   fprintf('OK\n');
catch err
   rethrow(err)
end
subplot(3,3,N)
show_slices(imgr);
eidors_msg('log_level',2);
eidors_cache on

function imdl = prepare_model
imdl = mk_common_model('b2C',16);



