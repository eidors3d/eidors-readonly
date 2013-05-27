% Reconstruct data on Gallery
% $Id$

n_iter=10;

gallery_3D_img.fwd_model.misc.compute_CCandSS='n';
for k= 1:n_iter
    eidors_msg('Iteration number %d',k,1);
    jacobian = calc_jacobian(gallery_3D_img);
    ref_data= fwd_solve(gallery_3D_img);
    residuals= real_data.meas-ref_data.meas;
    svj= svd(jacobian);
    % compute pseudo-inverse using only the largest singular values
    delta_params= pinv(jacobian,svj(1)/20.)*residuals;
    delta_params= delta_params.*gallery_3D_img.params_mapping.perturb;
    gallery_3D_img.params_mapping.params= gallery_3D_img.params_mapping.params + delta_params;
end

%% Solve final model and display results
ref_data= fwd_solve(gallery_3D_img);

subplot(211)
plot([ref_data.meas,real_data.meas]);
print_convert tutorial410c.png '-density 75';
