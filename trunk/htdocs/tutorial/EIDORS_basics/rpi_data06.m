% Absolute reconstructions $Id$

imdl = mk_common_model('b2c2',32);
imdl.fwd_model = fmdl;
imdl.reconst_type = 'absolute';
imdl.hyperparameter.value = 808;
imdl.solve = @GN_abs_solve;

for iter = [1,2,3];
   imdl.parameters.max_iteration = iter;
   img = inv_solve(imdl , vi);
   img.calc_colours.cb_shrink_move = [0.5,0.8,0];
   show_fem(img,[1,1]); axis off; axis image

   print_convert(sprintf('rpi_data06%c.png', 'a'-1+iter),'-density 60');
end
