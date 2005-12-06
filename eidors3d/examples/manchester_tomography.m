% Example to show reconstructions from
% the Process Tomography 
% group from U.Manchester, UK
%
% (C) 2005 by Stephen Murphy. Licensed under GPL version 2.
% $Id: manchester_tomography.m,v 1.1 2005-12-06 14:50:10 aadler Exp $

load 2planes_2rods_Opp
fwd_m1=eidors_obj('fwd_model',fwd_mdl);
fwd_m1=eidors_obj('fwd_model',fwd_mdl); 
fwd_m1.misc.perm_sym = '{y}';
img=eidors_obj('image','homog');


img.fwd_model= fwd_m1;
img.elem_data= mat_ref;
vh= fwd_solve(img);
%display_meas(fwd_m1)

%img.elem_data= ... 
%   drawsphere(fwd_m1.nodes,fwd_m1.elems,mat_ref,[0,0,.05],.05,.05);
%vi= fwd_solve(img);

if 0
fwd_m2= fwd_m1;
renumber=[1:2:16;18:2:32];renumber=renumber(:);
fwd_m2.electrode= fwd_m1.electrode(renumber);
fwd_m2.stimulation= mk_stim_patterns(16,1,'{ad}','{ad}',[],1);
img2=eidors_obj('image','homog');
img2.elem_data= mat_ref;
img2.fwd_model= fwd_m2;
vh2= fwd_solve(img2);
display_meas(fwd_m2)
end

fwd_m1.jacobian= 'np_calc_jacobian';
fwd_m1.

imdl= eidors_obj('inv_model','NP mdl');
imdl.solve= 'np_inv_solve';
%imdl.solve= 'aa_inv_conj_grad';
imdl.hyperparameter.value = 1e-3;
imdl.R_prior.func= 'np_calc_image_prior';
imdl.np_calc_image_prior.parameters= [3 1];
imdl.reconst_type= 'difference';
imdl.fwd_model= fwd_m1;

meas.type='data';
reference.type='data';
imr= inv_solve(imdl, meas, reference);
%imr= inv_solve(imdl, vi, vh);

img.elem_data= mat_ref;
J0 = np_calc_Jacobian(fwd_m1 , img);
img.elem_data= 100*mat_ref;
J1 = np_calc_Jacobian(fwd_m1 , img);

%hess = J'*J;
%[sol] = pcg(hess, J'*(meas.meas-reference.meas), 1e-4, 50);
%imr.elem_data=sol;
show_slices(imr, linspace(.01,.09,4)'*[inf,inf,1])

