% Compare different 2D reconstructions

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: compare_2d_algs.m,v 1.6 2005-12-05 17:55:59 aadler Exp $

imb=  mk_common_model('c2c',16);
e= size(imb.fwd_model.elems,1);
sigma= ones(e,1);
img= eidors_obj('image','');
img.elem_data= sigma;
img.fwd_model= imb.fwd_model;
vh= fwd_solve( img );

%sigma([65,81,82,101,102,122])=2; round
sigma([25,37,49:50,65:66,81:83,101:103,121:124])=2;
sigma([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=2;
img.elem_data= sigma;
vi= fwd_solve( img );

sig= sqrt(norm(vi.meas - vh.meas));
m= size(vi.meas,1);
vi.meas = vi.meas + .001*sig*randn(m,1);
figure(2); img.elem_data= img.elem_data - 1; show_slices(img); figure(1);

%show_slices(img);
imb=  mk_common_model('b2c',16);
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.fwd_model= imb.fwd_model;
inv2d.fwd_model.misc.perm_sym= '{y}';

switch 6
   case 1,
     inv2d.hyperparameter.value = 1e-3;
     inv2d.solve=       'aa_inv_solve';
     inv2d.image_prior.func= 'laplace_image_prior';

   case 2,
     inv2d.hyperparameter.value = 1e-3;
     inv2d.image_prior.func= 'laplace_image_prior';
     inv2d.solve=       'np_inv_solve';

   case 3,
     inv2d.hyperparameter.func = 'aa_calc_noise_figure';
     inv2d.hyperparameter.noise_figure= 2;
     inv2d.hyperparameter.tgt_elems= 1:4;
     inv2d.image_prior.func= 'aa_calc_image_prior';
     inv2d.solve=       'aa_inv_solve';

   case 4,
     inv2d.hyperparameter.value = [1e-2, 1e-6];
     inv2d.parameters.max_iterations= 5;
     inv2d.image_prior.func= 'ab_calc_tv_prior';
     inv2d.solve=       'ab_tv_diff_solve';

   case 5,
     inv2d.hyperparameter.value = 1e-2;
     inv2d.solve=       'aa_inv_total_var';
     inv2d.image_prior.func= 'laplace_image_prior';
     inv2d.parameters.max_iterations= 10;

   case 6,
     subplot(141); show_slices(img);
     inv2d.hyperparameter.value = 1e-4;
     inv2d.solve=       'aa_inv_total_var';
     inv2d.image_prior.func= 'laplace_image_prior';
     inv2d.parameters.max_iterations= 1;
     subplot(142); show_slices( inv_solve( inv2d, vi, vh) );
     inv2d.parameters.max_iterations= 3;
     subplot(143); show_slices( inv_solve( inv2d, vi, vh) );
     inv2d.parameters.max_iterations= 5;
     subplot(144); show_slices( inv_solve( inv2d, vi, vh) );
     return;

   otherwise,
     error('action unknown');
end

% 
% Step 3: Reconst and show image
% 
imgr= inv_solve( inv2d, vi, vh);
for i =1:length(imgr);
    ii= imgr(i).elem_data;
    ii=ii-1;
    ii=ii/max(ii(:));
    imgr(i).elem_data=ii;
end

figure(1); show_slices(imgr);
