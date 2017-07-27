%Make an ivnerse model - standard foward model inside
%imdl = mk_common_model('c2C',16);
%fmdl = imdl.fwd_model; %Extract model

for r=[0.5]
    sb = sprintf('solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;%0.3f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;',r);
    extra={'ball',sb};

    fmdl= ng_mk_cyl_models(0,[16],[0.2,0,0.05],extra); 
    ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.2).^2 + (ctr(:,2)-0.2).^2;
    img= mk_image(fmdl, 1 + 0.1*(ctr<0.2^2));

%    figure; show_fem(img);

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl = fix_model(fmdl);
fmdl.approx_type='tri3';

%Background, upper, lower conductivity and range
cond_bkg = 1; cond_lower = 0.2; cond_upper = 5;
cond_vals = linspace(cond_lower,cond_upper,100);

%Pixk group of pixels to perturb
%pixel_group = [102,123,103,83,66,82]; %b2C
%pixel_group = [327,364,328,292,259,291]; %c2C
pixel_group=fmdl.mat_idx{2}
%pixel_group = 1:length(fmdl.elems(:,1));

%Make an image on unit conductivity and solve
img = mk_image(fmdl,cond_bkg); img_b=img; 
vh = fwd_solve(img_b); %figure; plot(vh.meas); 
img.elem_data([pixel_group])=2;
%figure; show_fem(img,[1,0,1]);

%Calculate Jacobian and see difference between linear approximation
J=calc_jacobian(img_b);
J_pixel_group = J(:,pixel_group);
J_pixel_group_sum = sum(J(:,pixel_group),2);

%Phessian Jacobian
img.fwd_solve.get_all_meas=1;
d_k = fwd_solve(img); % doesn't contribute to cts as cached value here
u0= d_k.volt;
DU0 = calc_grad_potential(img, u0);
[~,~,~,J_phess] = calc_phessian_obj(fmdl,img,DU0,d_k.meas,'disc');
J_phess_pixel_group = J_phess(:,pixel_group);
J_phess_pixel_group_sum = sum(J_phess(:,pixel_group),2);

%Calculate Hessian and see difference between quadratic approximation
img.fwd_model.approx_type='tri3';
H = calc_hessian(img_b.fwd_model,img_b,pixel_group);
H_diag = calc_hessian_diag(img_b.fwd_model,img_b,pixel_group);


%Put into matrix
H_pixel_group = H(:,1:length(pixel_group),1:length(pixel_group));
H_pixel_group_diag = H_diag(:,1:length(pixel_group),1:length(pixel_group));

%figure; subplot(121); imagesc(squeeze(H_diag(1,:,:)));
%subplot(122); imagesc(squeeze(H(1,:,:)))

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
    dvlp{jj} = vh.meas + J_phess_pixel_group_sum*(cond_vals(jj)-cond_bkg);
    norm_vh_linp(jj)= norm(dvlp{jj},2);            
    dvq{jj} = vh.meas + J_pixel_group_sum*(cond_vals(jj)-cond_bkg) + ...
        0.5*H_pixel_group_sum*(cond_vals(jj)-cond_bkg)*(cond_vals(jj)-cond_bkg);
    norm_vh_quad(jj)= norm(dvq{jj},2);           
    
%     %% TODO ADD POLARISATION TENSOR GRADIENT/HESSIAN HERE
%     
%     [~, ~, C0] = calc_phessian_obj(img.fwd_model, img, DU0, delta_d, 'disc' );
%     H_phess = spdiags(C0,0,length(C0), length(C0)) + J_phess.'*J_phess;
%     H_phess_pixel_group_sum = sum(H_phess(:,pixel_group),2);    
%     dvqp{jj} = vh.meas + J_pixel_group_sum*(cond_vals(jj)-cond_bkg) + ...
%         0.5*H_phess_pixel_group_sum*(cond_vals(jj)-cond_bkg)*(cond_vals(jj)-cond_bkg);
    
end

%Plot the norm of voltage differences as function of contrast
figure; plot(cond_vals,norm_vi,'g*');
hold on; plot(cond_vals,norm_vh_lin,'r*');
hold on; plot(cond_vals,norm_vh_linp,'k*');
hold on; plot(cond_vals,norm_vh_quad,'b*');
legend('Non-linear','Linear approximation','Linear Phess','Quadratic approximation')
ylabel('|| V0 + DVh +D2Vh2 ||'); 
xlabel('Pertubation in ith pixel');

%Plot norm of voltage differences of linear and quadratic away from non-linear as function of contrast 
figure; plot(cond_vals,abs(norm_vi-norm_vh_lin)./abs(norm_vi),'r*');
hold on; plot(cond_vals,abs(norm_vi-norm_vh_linp)./abs(norm_vi),'k*');
hold on; plot(cond_vals,abs(norm_vi-norm_vh_quad)./abs(norm_vi),'b*');
legend('Linear-Nonlinear','Linear Phess - Nonlinear','Quadratic-Nonlinear')
ylabel('( ||Vi||-||Vh + DVh + D2Vh2||)'); 
xlabel('Pertubation in ith pixel');

end


%{
%Compute the Hessian ob objective
for jj=1:length(cond_vals)    
    fprintf(1,'Conductivity value %i of %i\n',jj,length(cond_vals));
    
    img.elem_data([pixel_group]) = cond_vals(jj);
    vi = fwd_solve(img); %Image with pixel perturbation        
    [H_obj,GN_only_obj,H_only_obj] = calc_hessian_obj(img_b.fwd_model,img_b,1:size(fmdl.elems,1),vi.meas-vh.meas);            
    H_only_obj_norm(jj) = norm(H_only_obj);
    GN_only_obj_norm(jj) = norm(GN_only_obj);
end

%Plot norm of voltage differences of linear and quadratic away from non-linear as function of contrast 
figure; plot(cond_vals,abs(GN_only_obj_norm-H_only_obj_norm)./abs(GN_only_obj_norm),'r*');
legend('||GNonly hess|| - ||Honly hess||')
ylabel('||GNonly hess|| - ||Honly hess||'); 
xlabel('Pertubation in ith pixel');
%}


