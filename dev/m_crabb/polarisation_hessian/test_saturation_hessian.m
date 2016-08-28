%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('c2C',16);
fmdl = imdl.fwd_model; %Extract model

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

%Background, upper, lower conductivity and range
cond_bkg = 1; cond_lower = 1; cond_upper = 3;
cond_vals = linspace(cond_lower,cond_upper,100);

%Pixk group of pixels to perturb
pixel_group = [327,364,328,292,259,291];
%pixel_group = [311,364];

%Make an image on unit conductivity and solve
img = mk_image(fmdl,cond_bkg); img_b=img; 
vh = fwd_solve(img_b); %figure; plot(vh.meas); 
img.elem_data([pixel_group])=2;
figure; show_fem(img,[1,0,1]);

%Calculate Jacobian and see difference between linear approximation
J=calc_jacobian(img_b);
J_pixel_group = J(:,pixel_group);
J_pixel_group_sum = sum(J(:,pixel_group),2);

%Calculate Hessian and see difference between quadratic approximation
img.fwd_model.approx_type='tri3';
H = calc_hessian(img_b.fwd_model,img_b,pixel_group);
H_pixel_group = H(:,1:length(pixel_group),1:length(pixel_group));
%for ff=1:size(H_pixel_group,1)
%    H_pixel_group_sum(ff,1) = ones(length(pixel_group),1)'*...
%        reshape(H_pixel_group(ff,:,:),length(pixel_group),length(pixel_group))*ones(length(pixel_group),1);
%end
H_pixel_group_sum = sum(sum( H_pixel_group ,2 ),3);

%Non-linear voltage difference
for jj=1:length(cond_vals)
    fprintf(1,'Conductivity value %i of %i\n',jj,length(cond_vals));
    img.elem_data([pixel_group]) = cond_vals(jj);
    vi = fwd_solve(img); %Image with pixel perturbation
    norm_vi(jj) = norm(vi.meas,2);
    dvl{jj} = vh.meas + J_pixel_group_sum*(cond_vals(jj)-cond_bkg);
    norm_vh_lin(jj)= norm(dvl{jj},2);        
    dvq{jj} = vh.meas + J_pixel_group_sum*(cond_vals(jj)-cond_bkg) + ...
        0.5*H_pixel_group_sum*(cond_vals(jj)-cond_bkg)*(cond_vals(jj)-cond_bkg);
    norm_vh_quad(jj)= norm(dvq{jj},2);             
end

%Plot the norm of voltage differences as function of contrast
figure; plot(cond_vals,norm_vi,'g*');
hold on; plot(cond_vals,norm_vh_lin,'r*');
hold on; plot(cond_vals,norm_vh_quad,'b*');
legend('Non-linear','Linear approximation','Quadratic approximation')
ylabel('|| Vi - V1||/||Vi||'); 
xlabel('Pertubation in ith pixel');

%Plot norm of voltage differences of linear and quadratic away from non-linear as function of contrast 
figure; plot(cond_vals,abs(norm_vi-norm_vh_lin)./abs(norm_vi),'r*');
hold on; plot(cond_vals,abs(norm_vi-norm_vh_quad)./abs(norm_vi),'b*');
legend('Linear-Nonlinear','Quadratic-Nonlinear')
ylabel('( ||Vi||-||Vh + DVh + D2Vh2||)'); 
xlabel('Pertubation in ith pixel');
