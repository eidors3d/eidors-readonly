function imgr= compare_3d_algs( algno )
% Compare different 3D reconstructions
% 
% algno=1     np_inv_solve            np_calc_image_prior
% algno=2     laplace_image_prior     np_inv_solve
% algno=3     ab_calc_tv_prior        np_inv_solve
% algno=4     ab_calc_tv_prior        ab_tv_diff_solve

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: compare_3d_algs.m,v 1.18 2007-08-29 09:14:09 aadler Exp $

calc_colours('ref_level','auto');

imb=  mk_common_model('n3r2',16);
e= size(imb.fwd_model.elems,1);
sigma= ones(e,1);
img= eidors_obj('image','');
img.elem_data= sigma;
img.fwd_model= imb.fwd_model;
vh= fwd_solve( img );

load datacom A B;
sigma(A)= 1.2; sigma(B)=0.8;
clear A B;
img.elem_data= sigma;
vi= fwd_solve( img );

sig= sqrt(norm(vi.meas - vh.meas));
m= size(vi.meas,1);
vi.meas = vi.meas + .001*sig*randn(m,1);

%show_slices(img);
inv3d= eidors_obj('inv_model', 'EIT inverse');
inv3d.reconst_type= 'difference';
inv3d.jacobian_bkgnd.value = 1;
inv3d.fwd_model= imb.fwd_model;
inv3d.fwd_model.misc.perm_sym= '{y}';

     iidx=1;
switch algno
   case 1,
     inv3d.hyperparameter.value = 1e-4;
     inv3d.solve=       'np_inv_solve';
     inv3d.R_prior=     'np_calc_image_prior';
     inv3d.np_calc_image_prior.parameters= [3 1]; %  deg=1, w=1

   case 2,
     inv3d.hyperparameter.value = 1e-3;
     inv3d.RtR_prior=    'laplace_image_prior';
     inv3d.solve=        'np_inv_solve';

   case 3,
     inv3d.hyperparameter.value = 1e-2;
     inv3d.R_prior=      'ab_calc_tv_prior';
     inv3d.solve=        'np_inv_solve';

   case 4,
     inv3d.hyperparameter.value = 1e-2;
     inv3d.ab_calc_tv_prior.alpha2 = 1e-5;
     inv3d.parameters.max_iterations= 5;
     inv3d.parameters.term_tolerance= 1e-3;
     inv3d.parameters.keep_iterations= 1;
     inv3d.R_prior=      'ab_calc_tv_prior';
     inv3d.solve=        'ab_tv_diff_solve';
     iidx=[1,2,5];

   case 5,
     inv3d.hyperparameter.value = 1e-3;
     inv3d.R_prior=      'ab_calc_tv_prior';
     inv3d.solve=        'aa_inv_total_var';
     inv3d.parameters.max_iterations= 5;

   otherwise,
     error('action unknown');
end

% 
% Step 3: Reconst and show image
% 
imgr= inv_solve( inv3d, vi, vh);
show_slices(imgr(iidx), [.5:1:2.5]'*[Inf,Inf,1]);

