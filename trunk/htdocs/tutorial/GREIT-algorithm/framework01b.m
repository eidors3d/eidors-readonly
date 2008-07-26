% fwd_model $Id: framework01.m 1536 2008-07-26 15:53:40Z aadler $

% CALCULATE MODEL CORRESPONDENCES
[rmdl,c2f] = mk_grid_model(ng_mdl_16x1_fine, xvec, yvec, zvec);

% CALCULATE JACOBIAN AND SAVE IT

img= eidors_obj('image','GREIT-ng_mdl');
img.fwd_model= ng_mdl_16x1_fine;
img.fwd_model.coarse2fine = c2f;
img.rec_model= rmdl;
img.elem_data= ones(size(img.fwd_model,1));

% ADJACENT STIMULATION PATTERNS
img.fwd_model.stimulation= mk_stim_patterns(16, 1, ...
             [0,1],[0,1], {'do_redundant', 'no_meas_current'}, 1);

% SOLVERS
img.fwd_model.system_mat= @aa_calc_system_mat;
img.fwd_model.solve=      @aa_fwd_solve;
img.fwd_model.jacobian=   @aa_calc_jacobian;

J= calc_jacobian(img);
save GREIT_Jacobian_ng_mdl_fine J
