% Compare different 2D reconstructions

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: compare_2d_algs.m,v 1.2 2005-12-01 18:20:10 aadler Exp $

imb=  mk_common_model('b2c',16);
e= size(imb.fwd_model.elems,1);
sigma= ones(e,1);
img= eidors_obj('image','');
img.elem_data= sigma;
img.fwd_model= imb.fwd_model;
vh= fwd_solve( img );

%sigma([65,81,82,101,102,122])=2; round
sigma([25,37,49:50,65:66,81:83,101:103,121:124])=2;
img.elem_data= sigma;
vi= fwd_solve( img );

sig= sqrt(norm(vi.meas - vh.meas));
m= size(vi.meas,1);
vi.meas = vi.meas + .000001*sig*randn(m,1);

%show_slices(img);
inv2d.name= 'EIT inverse';
%inv2d.solve=       'aa_inv_solve';
 inv2d.solve=       'np_inv_solve'; iidx=1;
%inv2d.solve=       'ab_tv_diff_solve'; iidx= [2:10];
%inv2d.solve=       'aa_inv_total_var';
 inv2d.hyperparameter.value = 1e-3;
%inv2d.hyperparameter.func = 'aa_calc_noise_figure';
%inv2d.hyperparameter.noise_figure= 2;
%inv2d.hyperparameter.tgt_elems= 1:4;
%inv2d.image_prior.func= 'laplace_image_prior';
 inv2d.image_prior.func= 'ab_calc_tv_prior';
%inv2d.image_prior.func= 'aa_calc_image_prior';
inv2d.reconst_type= 'difference';
inv2d.fwd_model= imb.fwd_model;
inv2d.fwd_model.misc.perm_sym= '{y}';
inv2d= eidors_obj('inv_model', inv2d);

% 
% Step 3: Reconst and show image
% 
img= inv_solve( inv2d, vi, vh);
show_slices(img(iidx));

