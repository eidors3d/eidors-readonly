% Reconstruct data on Gallery
% $Id: tutorial410c.m,v 1.1 2007-06-15 19:56:56 aadler Exp $

n_iter=9;

gallery_3D_img.fwd_model.misc.compute_CCandSS='n';
for k= 1:n_iter
    eidors_msg('Iteration number %d',k,1);
    jacobian = dg_calc_jacobian(gallery_3D_img);
    [ref_data,gallery_3D_img]= dg_fwd_solve(gallery_3D_img);
    residuals= real_data.meas-ref_data.meas;
    disp('Solving the inverse problem');
    svj= svd(jacobian);
    % compute pseudo-inverse using only the largest singular values
    delta_params= pinv(jacobian,svj(1)/20.)*residuals;
    delta_params= delta_params.*gallery_3D_img.params_mapping.perturb;
    gallery_3D_img.params_mapping.params= gallery_3D_img.params_mapping.params + delta_params;
end

%% Solve final model and display results
[ref_data,gallery_3D_img]= dg_fwd_solve(gallery_3D_img);

subplot(211)
plot([ref_data.meas,real_data.meas]);
 print -r75 -dpng tutorial410c.png;
