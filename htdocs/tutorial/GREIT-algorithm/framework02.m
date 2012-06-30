% fwd_model $Id$

% CALCULATE JACOBIAN AND SAVE IT

img= mk_image(fmdl, 1);
img.fwd_model.coarse2fine = c2f;
img.rec_model= rmdl;

% ADJACENT STIMULATION PATTERNS
img.fwd_model.stimulation= mk_stim_patterns(16, 1, ...
             [0,1],[0,1], {'do_redundant', 'no_meas_current'}, 1);

% SOLVERS
img.fwd_model.system_mat= @system_mat_1st_order;
img.fwd_model.solve=      @fwd_solve_1st_order;
img.fwd_model.jacobian=   @jacobian_adjoint;

J= calc_jacobian(img);

map = reshape(sum(c2f,1),pixel_grid,pixel_grid)>0;
save GREIT_Jacobian_ng_mdl_fine J map
