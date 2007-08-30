% Solve resistor model
% $Id: tutorial020c.m,v 1.2 2007-08-30 03:58:28 aadler Exp $

% Now we complete the fwd_model
r_mdl.jacobian= @perturb_jacobian;

% Now create an inverse model
i_mdl= eidors_obj('inv_model','resistor inverse');
i_mdl.fwd_model= r_mdl;
i_mdl.jacobian_bkgnd.value= 1000;

% regulatization not needed for this problem
i_mdl.RtR_prior= @tikhonov_image_prior;
i_mdl.hyperparameter.value= 0;

i_mdl.reconst_type= 'difference';
i_mdl.solve= @aa_inv_solve;

% Reconstruct resistor change
reconst= inv_solve(i_mdl, data_1k0, data_1k2);
