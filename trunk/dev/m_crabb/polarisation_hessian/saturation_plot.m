%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('c2C',16);
fmdl = imdl.fwd_model; %Extract model

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

%Make an image and show image
img = mk_image(fmdl,1);
img_b=img;
figure; show_fem(img,[0,1,3])

%Solve forward model
vh = fwd_solve(img);
figure; plot(vh.meas); 

%Pixk group of pixels to perturb
pixel_group = [327,364,328,292,259,291];

%Conductivity values
cond_vals = logspace(0,.5,10);

for jj=1:length(cond_vals)
    jj;
    img.elem_data([pixel_group]) =cond_vals(jj);
    vi = fwd_solve(img); %Image with pixel perturbation
    norm_vi_diff(jj) = norm(vi.meas-vh.meas);
    
end


%Calculate Jacobian and see difference between linear approximation
J=calc_jacobian(img_b);
J_pixel_group = J(:,pixel_group);
J_pixel_group_sum = sum(J(:,pixel_group),2);

for jj=1:length(cond_vals)
    jj;
    dv =J_pixel_group_sum*(cond_vals(jj)-1);
    norm_vi_diff_lin(jj)= norm(dv);                
end

%Calculate Hessian and see difference between quadratic approximation
img.fwd_model.approx_type='tri3';
H = calc_hessian_select(img_b.fwd_model,img_b,pixel_group);

H_pixel_group = H(:,1:length(pixel_group),1:length(pixel_group));
H_pixel_group_sum = sum(sum( H_pixel_group ,2 ),3);

for jj=1:length(cond_vals)
    jj;
    dv =J_pixel_group_sum*(cond_vals(jj)-1) + 0.5*H_pixel_group_sum*(cond_vals(jj)-1)*(cond_vals(jj)-1);
    norm_vi_diff_quad(jj)= norm(dv);                
end


figure; semilogx(cond_vals,norm_vi_diff,'r*');
hold on; semilogx(cond_vals,norm_vi_diff_lin,'b*');
hold on; semilogx(cond_vals,norm_vi_diff_quad,'g*');
legend('Non-linear perturbation','Linear approximation','Quadratic approximation')
ylabel('|| Vi - V0||'); 
xlabel('Pertubation in ith pixel');

