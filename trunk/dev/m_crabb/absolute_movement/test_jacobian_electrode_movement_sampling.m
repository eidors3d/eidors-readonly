clc; close all;
run ~/EIT/Code/mk_paths.m

n_dim=2; 

if(n_dim==2)
%No .of electrodes and the anomaly parameter
nelec=16; elec_t=[0.05,0,0.025]; %CEM

%Homogeneous model
%Create circle homogeneous CEM (0.05) elecs equal spaced
elec_pos=zeros(nelec,2); elec_pos_cart=zeros(nelec,2);
for i=1:nelec
    elec_pos(i,1)=(i-1)*360/nelec; elec_pos(i,2)=0;
    elec_pos_cart(i,1)=-sin(elec_pos(i,1)*pi/180);
    elec_pos_cart(i,2)=cos(elec_pos(i,1)*pi/180);
end

mdl_h=ng_mk_cyl_models([0,1,0.1],elec_pos,elec_t);

%Make some stimulation patterns and add to models
stim=mk_stim_patterns(16,1,[0 1],[0 1]);
mdl_h.stimulation= stim; 

%Make inhomogeneous and homogeneous image
img_h=mk_image(mdl_h,1);
figure; show_fem(img_h,[1,0,0]);

%Homogeneous and inhomgeneous voltages
v_h=fwd_solve(img_h.fwd_model,img_h);
elec_comp = calc_electrode_components(mdl_h);
%Get the Movement Jacobian and the components into normal and tangential
tic; [Jold] = jacobian_movement_eidors_electrode_tangential(img_h.fwd_model,img_h); toc;
tic; [Jnew] = jacobian_movement_sampling_electrode_tangential(img_h.fwd_model,img_h); toc;

elseif(n_dim==3)
%No .of electrodes and the anomaly parameter
nelec=16; elec_t=[0.05,0,0.05]; %CEM

%Homogeneous model
%Create circle homogeneous CEM (0.05) elecs equal spaced
elec_pos=zeros(nelec,2); elec_pos_cart=zeros(nelec,2);
for i=1:nelec
    elec_pos(i,1)=(i-1)*360/nelec; elec_pos(i,2)=0.5;
end
mdl_h=ng_mk_cyl_models([1,1,0.1],elec_pos,elec_t);

%Make some stimulation patterns and add to models
stim=mk_stim_patterns(16,1,[0 1],[0 1]);
mdl_h.stimulation= stim; 
img_h=mk_image(mdl_h,1);
figure; show_fem(img_h,[1,0,0]);

%Get the electrode components
elec_comp = calc_electrode_components(mdl_h);

%Get the Movement Jacobian and the components into normal and tangential
[Jold] = jacobian_movement_eidors_electrode_tangential(img_h.fwd_model,img_h);
%[Jnew] = jacobian_movement_sampling_electrode_tangential(img_h);

end
