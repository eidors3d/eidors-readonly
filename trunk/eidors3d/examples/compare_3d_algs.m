% Compare different 3D reconstructions

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: compare_3d_algs.m,v 1.3 2005-12-02 11:49:57 aadler Exp $

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
vi.meas = vi.meas + .01*sig*randn(m,1);

%show_slices(img);
inv3d= eidors_obj('inv_model', 'EIT inverse');
inv3d.reconst_type= 'difference';
inv3d.fwd_model= imb.fwd_model;
inv3d.fwd_model.misc.perm_sym= '{y}';

     iidx=1;
switch 3
   case 1,
     inv3d.hyperparameter.value = 1e-4;
     inv3d.solve=            'np_inv_solve';
     inv3d.image_prior.func= 'np_calc_image_prior';
     inv3d.image_prior.parameters= [3 1]; %  deg=1, w=1

   case 2,
     inv3d.hyperparameter.value = 1e-4;
     inv3d.image_prior.func=  'laplace_image_prior';
     inv3d.solve=             'np_inv_solve';

   case 3,
     inv3d.hyperparameter.value = 1e-2;
     inv3d.image_prior.func=  'ab_calc_tv_prior';
     inv3d.solve=             'np_inv_solve';

   otherwise,
     error('action unknown');
end

% 
% Step 3: Reconst and show image
% 
imgr= inv_solve( inv3d, vi, vh);
show_slices(imgr(iidx), [.5:1:2.5]'*[Inf,Inf,1]);

