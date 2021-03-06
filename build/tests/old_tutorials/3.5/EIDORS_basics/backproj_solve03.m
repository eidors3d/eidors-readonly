% $Id$

tutorial120a; % get the model from a related tutorial

% Gauss Newton Solver
inv_GN= eidors_obj('inv_model','GN_solver','fwd_model', img.fwd_model);
inv_GN.reconst_type= 'difference';
inv_GN.solve= @np_inv_solve;
inv_GN.RtR_prior= @noser_image_prior;
inv_GN.jacobian_bkgnd.value= 1;
inv_GN.hyperparameter.value= 0.03;

imgr= inv_solve(inv_GN, vh,vi);
imgr.calc_colours.ref_level=0;
subplot(131); show_fem(imgr);
axis equal; axis off

% Backprojection Solver
inv_BP= eidors_obj('inv_model','BP_solver','fwd_model', img.fwd_model);
inv_BP.reconst_type= 'difference';
inv_BP.solve= @backproj_solve;
inv_BP.backproj_solve.type= 'naive';

imgr= inv_solve(inv_BP, vh,vi);
imgr.calc_colours.ref_level=0;
subplot(132); show_fem(imgr);
axis equal; axis off

inv_BP.backproj_solve.type= 'simple_filter';

imgr= inv_solve(inv_BP, vh,vi);
imgr.calc_colours.ref_level=0;
subplot(133); show_fem(imgr);
axis equal; axis off

print_convert backproj_solve03a.png
