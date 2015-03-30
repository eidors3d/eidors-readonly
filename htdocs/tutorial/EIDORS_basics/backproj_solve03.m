% $Id$

tutorial120a; % get the model from a related tutorial

% Gauss Newton Solver
inv_GN= eidors_obj('inv_model','GN_solver','fwd_model', img.fwd_model);
inv_GN.reconst_type= 'difference';
inv_GN.solve= @inv_solve_diff_GN_one_step;
inv_GN.RtR_prior= @prior_noser;
inv_GN.jacobian_bkgnd.value= 1;
inv_GN.hyperparameter.value= 0.2;

imgr= inv_solve(inv_GN, vh,vi);
imgr.calc_colours.ref_level=0;
subplot(131); show_fem(imgr);
axis equal; axis off

% Backprojection Solver
inv_BP= eidors_obj('inv_model','BP_solver','fwd_model', img.fwd_model);
inv_BP.reconst_type= 'difference';
inv_BP.solve= @inv_solve_backproj;
inv_BP.inv_solve_backproj.type= 'naive';

imgr= inv_solve(inv_BP, vh,vi);
imgr.calc_colours.ref_level=0;
subplot(132); show_fem(imgr);
axis equal; axis off

inv_BP.inv_solve_backproj.type= 'simple_filter';

imgr= inv_solve(inv_BP, vh,vi);
imgr.calc_colours.ref_level=0;
subplot(133); show_fem(imgr);
axis equal; axis off

print_convert inv_solve_backproj03a.png
