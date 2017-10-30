fmdl= ng_mk_cyl_models([0,1,0.2],[8],[0.2,0,1]); 
show_fem(fmdl);
    
stim = mk_stim_patterns(8,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model
fmdl = fix_model(fmdl);
fmdl.approx_type='tri3';
nelems = size(fmdl.elems,1);

%uniform sigma \in [0.5,1.5]
img = mk_image(fmdl,0.5 + rand(size(fmdl.elems,1),1));
%figure; show_fem(img)




J=jacobian_adjoint(img.fwd_model,img)
Jp = jacobian_perturb(img.fwd_model,img)

tic; Hp = hessian_perturb(img.fwd_model,img); toc;
tic; H = calc_hessian(img.fwd_model,img,[1:nelems]); toc;

h_norm=0;h_diff_norm=0;
j_norm=0;j_diff_norm=0;
for i=1:size(J,1)
h_diff_norm = h_diff_norm + ( norm(squeeze(H(i,:,:)-Hp(i,:,:))) )^2;
h_norm = h_norm + ( norm(squeeze(H(i,:,:))) )^2;

j_diff_norm = j_diff_norm + ( norm(squeeze(J(i,:)-Jp(i,:))) )^2;
j_norm = j_norm + ( norm(squeeze(J(i,:))) )^2;
end


j_diff_norm/j_norm
j_norm

h_diff_norm/h_norm
h_norm