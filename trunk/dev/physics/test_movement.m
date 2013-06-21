function test_movement(test_no)
if nargin < 1
   test_no =  1:10;
end


mdlc = mk_common_model('c2c');
img = mk_image( mdlc, 1,'conductivity');
img.conductivity.elem_data([75,93,94,113,114,136]) = 1.2;
img.conductivity.elem_data([105,125,126,149,150,174]) = 0.8;
img.movement.electrode_data = zeros(length(img.fwd_model.electrode)*2,1);
img.physics_data_mapper = 'conductivity';
img.fwd_model.solve = @fwd_solve_elec_move;

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
         
            
      otherwise 
         fprintf('  (no more tests)\n')
         break
   end
end