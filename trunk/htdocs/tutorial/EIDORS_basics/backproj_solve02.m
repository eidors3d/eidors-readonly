% $Id: backproj_solve02.m,v 1.1 2007-10-15 18:09:43 aadler Exp $

inv_GN= eidors_obj('inv_model','GN_solver','fwd_model', img.fwd_model);
inv_GN.reconst_type= 'difference';
inv_GN.solve= @np_inv_solve;
inv_GN.RtR_prior= @noser_image_prior;
inv_GN.jacobian_bkgnd.value= 1;
inv_GN.hyperparameter.value= 1e-3;

imgr= inv_solve(inv_GN, vh,vi);
subplot(121)
show_fem(imgr);
axis equal; axis off



