%Reconstruction model
[cmdl]= mk_grid_model([], 2.5+[-50,-20,0:10:320,340,370], ...
                             -[0:2.5:10, 15:5:25,30:10:80,100,120]);
c2f = mk_coarse_fine_mapping( fmdl, cmdl);

fmdl.coarse2fine = c2f;
imdl= eidors_obj('inv_model','test');
imdl.fwd_model= fmdl;
imdl.rec_model= cmdl;
   imdl.reconst_type = 'difference';
   imdl.RtR_prior = @prior_laplace;
   imdl.solve = @inv_solve_diff_GN_one_step;
   imdl.hyperparameter.value = 0.3;
   imdl.fwd_model.normalize_measurements = 1;
   imdl.jacobian_bkgnd.value = 0.03;

hold on
show_fem(cmdl);
view([0 -0.2 1])
hold off
print_convert pont_pean01a.png
