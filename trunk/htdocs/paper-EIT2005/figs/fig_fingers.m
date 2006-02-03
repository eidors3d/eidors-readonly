% Generate fingers images for EIDORS paper
% $Id: fig_fingers.m,v 1.1 2006-02-03 18:34:07 aadler Exp $
imdl=mk_common_model('n3r2');
im=eidors_obj('image','3 fingers','fwd_model',imdl.fwd_model);
rr= ones(size(imdl.fwd_model.elems,1),1);
im.elem_data=rr;
vh= fwd_solve(im);

fingerelems=  [ ...
 232 233 234 238 239 240 253 254 255 ...
 508 509 510 514 515 516 529 530 531];
rr(fingerelems)=1.5;
im.elem_data=rr;
vi= fwd_solve(im);
calc_colours('ref_level',1);
show_fem(im)

imdl.hyperparameter.value= 1e-4;
irec= inv_solve(imdl,vi,vh);
calc_colours('ref_level',0);
show_fem(irec)
show_slices(irec,4)

