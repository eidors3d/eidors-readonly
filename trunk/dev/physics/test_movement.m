function test_movement(test_no)
if nargin < 1
   test_no =  1:20;
end


mdlc = mk_common_model('c2c');
img = mk_image( mdlc, 1,'conductivity');
vh = fwd_solve(img);
img.conductivity.elem_data([75,93,94,113,114,136]) = 1.2;
img.conductivity.elem_data([105,125,126,149,150,174]) = 0.8;
vi = fwd_solve(img);
img.movement.electrode_data = zeros(length(img.fwd_model.electrode)*2,1);
img.data_mapper = 'conductivity';
img.fwd_model.solve = @fwd_solve_elec_move;

mdl2dim = mk_common_model('b2c');
mdl2dim.hyperparameter.value= 0.01;

for i = test_no
   fprintf('Test %d\t:',i);
   switch i
      case 1
         data_mv = fwd_solve(img);
         data    = fwd_solve(rmfield(img,'movement'));
         
         assert(all(data_mv.meas == data.meas));
         fprintf('  OK\n');
      case 2
         data_mv = fwd_solve(img);
         img2    = rmfield(img,'movement');
         img2.fwd_model.solve = 'eidors_default';
         data    = fwd_solve(img2);
         
         assert(all(data_mv.meas == data.meas));
         fprintf('  OK\n');
         
              
      case 3
         subplot(121)
         show_fem(img,[0 1]);
         img.movement.electrode_data(5) = 0.1; % x coord of 5th elec
         data_mv = fwd_solve(img);
         
         img2    = rmfield(img,'movement');
         img2.fwd_model.solve = 'eidors_default';
         
         img2.fwd_model.nodes(img2.fwd_model.electrode(5).nodes,1) = ... 
            img2.fwd_model.nodes(img2.fwd_model.electrode(5).nodes,1) + 0.1;
         
         subplot(122)
         show_fem(img2,[0,1]);
         
         data    = fwd_solve(img2);
         assert(all(data_mv.meas == data.meas));
         fprintf('  OK\n');
         
      case 4
         
         subplot(121)
         show_fem(inv_solve(mdl2dim,vh,vi));
         
         subplot(122)
         img.movement.electrode_data(5) = 0.1; % x coord of 5th elec
         data_mv = fwd_solve(img);
         show_fem(inv_solve(mdl2dim,vh,data_mv));
         fprintf('  OK\n');
      case 5 
         img.fwd_model.jacobian = @jacobian_movement;
         J = calc_jacobian(img);
         e_move = zeros(32,1);
         e_move(5) = 0.1;
         params = [zeros(size(img.conductivity.elem_data)); e_move]; 
         % difference in data caused by movement according to Jacobian
         diff_J = J*params;
         
         % data before movement
         data = fwd_solve(img);
         img.movement.electrode_data(5) = 0.1; % x coord of 5th elec
         % data after movement
         data_mv = fwd_solve(img);
         
         diff_F = data_mv.meas - data.meas;
         
         clf
         plot(diff_F);
         hold all
         plot(diff_J);
         plot((diff_J-diff_F),'r','linewidth',2)
         legend('fwd\_solve diff', 'Jacobian diff', 'error');
         hold off
         
         fprintf('  OK\n');
      case 6
         mdlc = mk_common_model('c2c');
         img = mk_image( mdlc, 1,'conductivity');
         J = calc_jacobian(img);
         vh = fwd_solve(img);
         img.conductivity.elem_data([75,93,94,113,114,136]) = 1.2;
         img.conductivity.elem_data([105,125,126,149,150,174]) = 0.8;
         vi = fwd_solve(img);
         
         diff_J = J*(img.conductivity.elem_data-1);
         diff_F = vi.meas - vh.meas;
         
         clf
         plot(diff_F);
         hold all
         plot(diff_J);
         plot((diff_J-diff_F),'r','linewidth',2)
         legend('fwd\_solve diff', 'Jacobian diff', 'error');
         hold off
         
         fprintf('  OK\n');
      case 7
         mdlc.jacobian_bkgnd = img;
         R1 = calc_RtR_prior(mdlc);
         mdlc.RtR_prior = @prior_movement;
         R2 = calc_RtR_prior(mdlc);
         sz = size(R1);
         assert(all(all(R2(1:sz(1),1:sz(2)) == R1)));
         fprintf('  OK\n');
         
      case 8
         img = mk_image( mdlc, 1,'conductivity');
         img.movement.electrode_data = zeros(length(img.fwd_model.electrode)*2,1);
         img.data_mapper = 'conductivity';
         img.fwd_model.solve = @fwd_solve_elec_move;
         img.params_mapper = {'conductivity.elem_data', 'movement.electrode_data'};
         mdlc.jacobian_bkgnd = rmfield(img,'fwd_model');
         mdlc.fwd_model.jacobian = @jacobian_movement;
         mdlc.RtR_prior =          @prior_movement;
         mdlc.prior_movement.parameters = sqrt(1e2/1);
         
         % find the right scale
         mdlc.inv_solve_diff_GN_one_step.calc_step_size = 1;
         
         % Solve inverse problem for mdlM eidors_obj model.
         imgc = inv_solve(mdlc, vh, vi);
         fprintf('  OK\n');
         
      case 9
         img = mk_image( mdlc, 1,'conductivity');
         img.movement.electrode_data = zeros(length(img.fwd_model.electrode)*2,1);
         img.data_mapper = 'conductivity';
         img.fwd_model.solve = @fwd_solve_elec_move;
         img.params_mapper = {'conductivity.elem_data', 'movement.electrode_data'};
         mdlc.jacobian_bkgnd = rmfield(img,'fwd_model');
         mdlc.fwd_model.jacobian = @jacobian_movement;
         mdlc.RtR_prior =          @prior_movement;
         mdlc.prior_movement.parameters = sqrt(1e2/1);
         
         vh = fwd_solve(img);
         
         img.conductivity.elem_data([75,93,94,113,114,136]) = 1.2;
         img.conductivity.elem_data([105,125,126,149,150,174]) = 0.8;
         img.movement.electrode_data(5) = 0.1; % x coord of 5th elec
         vi = fwd_solve(img);         
         % find the right scale
         mdlc.inv_solve_diff_GN_one_step.step_size = 1;

         imgc = inv_solve(mdlc, vh, vi);
         show_fem(imgc);
         fprintf('  OK\n');
      case 10
         % absolute movement
         
         mdlc = mk_common_model('c2c');
         
         img = mk_image( mdlc, 1,'conductivity');
         img.movement.electrode_data = zeros(length(img.fwd_model.electrode)*2,1);
         img.data_mapper = 'conductivity';
         img.fwd_model.solve = @fwd_solve_elec_move;
         img.params_mapper = {'conductivity.elem_data', 'movement.electrode_data'};
         
         mdlc.jacobian_bkgnd = rmfield(img,'fwd_model');
         mdlc.fwd_model.jacobian = @jacobian_movement;
         mdlc.RtR_prior =          @prior_movement;
         mdlc.prior_movement.parameters = sqrt(1e2/1);
         
         img.conductivity.elem_data([75,93,94,113,114,136]) = 1.2;
         img.conductivity.elem_data([105,125,126,149,150,174]) = 0.8;
         img.movement.electrode_data(5) = 0.1; % x coord of 5th elec
         vi = fwd_solve(img);         
         
         mdlc.solve = 'inv_solve_abs_GN';
         mdlc.reconst_type = 'absolute';
         mdlc.inv_solve.max_iterations = 10;
         imgc = inv_solve(mdlc, vi);
         
         f = clf;
         subplot(131)
         show_fem_move(img);
         subplot(132)
         show_fem_move(imgc);
         
         mdlc.solve = 'inv_solve_abs_gn';
         mdlc.inv_solve.verbose = 0; % don't spit out figures
         imgc = inv_solve(mdlc, vi);
         
         subplot(133)
         show_fem_move(imgc);
         
         fprintf('  OK\n');
         
      case 11
         % old tutorials
         path = eidors_obj('eidors_path');
         path = [path '/../htdocs/tutorial/adv_image_reconst/'];
         tuts = {'move_2d01', 'move_2d02', 'move_3d01', 'move_3d02'};
         clf
         for i = 1:length(tuts)
            run([path tuts{i}]);
         end
      case 12
         % new tutorials
         path = [];
         tuts = {'move_2d01', 'move_2d02'};%, 'move_3d01', 'move_3d02'};
         clf
         for i = 1:length(tuts)
            run([path tuts{i}]);
         end

      otherwise 
         fprintf('  (no more tests)\n')
         break
   end
end
