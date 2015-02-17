%Choose shape, noise and movement and hyperparams
noise=50; shape='sphere'; move_fac=0.1;          
hpmt_params_diff=0.01; hp_params=[10^-5]; hpmt_params=[1250]; 
shape_iteration=8; max_its=5;
   
%Parameter for conductivity, electrode shapes
cond_inc=3; elec_ref=0.05; dom_ref=0.5; dom_h=1;
n_elec_m=14; n_elec_mt=9; h_mt=0.5; elec_t=[0.05,0,elec_ref]; 
n_elec_t=n_elec_m+2*n_elec_mt+2; n_elec=n_elec_t;

%Stimulation patterns
stim=mk_stim_patterns(n_elec_t,1,[0 1],[0 1]);

%% CREATE COARSE, FINE AND ACTUAL MODEL
elec_pos=zeros(n_elec_m+2*n_elec_mt+2,6); elec_pos_cart=elec_pos;
for i=1
    elec_pos(i,1)=0.0;
    elec_pos(i,2)=0.0;
    elec_pos(i,3)=1.0;
    elec_pos(i,4:6)=elec_pos(i,1:3);       
end
for i=1+1:n_elec_mt+1
    psi=acos(h_mt);     
    elec_pos(i,1)=cos((i-2)*2*pi/n_elec_mt)*sin(psi);
    elec_pos(i,2)=sin((i-2)*2*pi/n_elec_mt)*sin(psi);
    elec_pos(i,3)=h_mt;
    elec_pos(i,4:6)=elec_pos(i,1:3);          
end
for i=1+n_elec_mt+1:1+n_elec_mt+n_elec_m
    elec_pos(i,1)=cos((i-2-n_elec_mt)*2*pi/n_elec_m);
    elec_pos(i,2)=sin((i-2-n_elec_mt)*2*pi/n_elec_m);
    elec_pos(i,3)=0;
    elec_pos(i,4:6)=elec_pos(i,1:3);          
end
for i=1+n_elec_mt+n_elec_m+1:1+n_elec_mt+n_elec_m+n_elec_mt
    psi=acos(-h_mt);     
    elec_pos(i,1)=cos((i-2-n_elec_mt-n_elec_m)*2*pi/n_elec_mt)*sin(psi);
    elec_pos(i,2)=sin((i-2-n_elec_mt-n_elec_m)*2*pi/n_elec_mt)*sin(psi);
    elec_pos(i,3)=-h_mt;
    elec_pos(i,4:6)=elec_pos(i,1:3);          
end
for i=2+n_elec_m+2*n_elec_mt
    elec_pos(i,1)=0;
    elec_pos(i,2)=0;
    elec_pos(i,3)=-1.0;
    elec_pos(i,4:6)=elec_pos(i,1:3);          
end

%Homogeneous sphere - 1 top/bottom, 8 above/below, 16 middle
shape_str = ['solid top     = ellipsoid(0,0,0; 0,0,1; 1,0,0; 0,1,0); \n' ...
    'solid mainobj = top -maxh=0.8;\n'];
elec_shape=elec_t;
elec_obj = 'top';
mdl_h = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
figure; show_fem(mdl_h);
nodes=mdl_h.nodes;

%Homogeneous sphere with no electrodes 
shape_str_c = ['solid top     = ellipsoid(0,0,0; 0,0,1.05; 1.05,0,0; 0,1.05,0); \n' ...
    'solid mainobj = top -maxh=0.25;\n'];
mdl_c = ng_mk_gen_models(shape_str_c, [], [] , {} );
figure; show_fem(mdl_c);
nodes_c=mdl_c.nodes;

%Move the electrode positions and renormalize
elec_pos_i=elec_pos;
for i=1:length(elec_pos_i(:,1))
    for j=1:3
        elec_pos_i(i,j)= elec_pos(i,j) + move_fac*(-1+2*rand(1,1)); 
    end
    elec_pos_i(i,1:3)=elec_pos_i(i,1:3)/norm(elec_pos_i(i,1:3));
    elec_pos_i(i,4:6)=elec_pos_i(i,1:3);
end
%Inhomogeneous sphere - 1 top/bottom, 8 above/below, 16 middle
shape_str_i = ['solid top     = ellipsoid(0,0,0; 0,0,1; 1,0,0; 0,1,0); \n' ...
               'solid ball    = sphere(0.0,0.3,0.3;0.3); tlo ball; \n' ...
               'solid mainobj = top and not ball -maxh=0.8; \n'];
elec_shape_i=elec_t;
elec_obj_i = 'top';
mdl_i = ng_mk_gen_models(shape_str_i, elec_pos_i, elec_shape_i, elec_obj_i);
figure; show_fem(mdl_i);
nodes=mdl_i.nodes;

%Make inhomogeneous image
img_h=mk_image(mdl_h,1);
img_i=mk_image(mdl_i,1); img_i.elem_data(mdl_i.mat_idx{1})=cond_inc; 
mdl_h.stimulation= stim; img_h.fwd_model.stimulation=stim;
mdl_i.stimulation= stim; img_i.fwd_model.stimulation=stim;

%Inhomgeneous voltages
v_i=fwd_solve(img_i.fwd_model,img_i);
v_h=fwd_solve(img_h.fwd_model,img_h); v_hd = v_h.meas; 
v_in=add_noise(noise,v_i); v_i=v_in; %Add noise

%% STEP 1 : ELECTRODES ONLY
for kk=1:shape_iteration
%Create the inverse model
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
mdl_h.jacobian=@jacobian_adjoint; 
img_h.fwd_model.jacobian=@jacobian_adjoint;
inv_tik_2d.fwd_model= mdl_h;
inv_tik_2d.reconst_type='absolute';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.RtR_prior=@prior_movement_tangential_only; 
inv_tik_2d.solve=@inv_solve_diff_tangential;
inv_tik_2d.hyperparameter.value(1)=hpmt_params_diff;   
inv_tik_2d.parameters.max_iterations=1;

%Add the  information on the type
inv_tik_2d.shape=shape;

%Add some new fields to the inversion
inv_tik_2d.img_fine=img_h;
inv_tik_2d.mdl_coarse=mdl_c;

%Differene solve on the inverse model (this has mdl_fine and c2f)
img_recon=inv_solve(inv_tik_2d,v_i);
%Calculate the new electrode positions from tangent data
elec_comp_h=calc_electrode_components(img_h.fwd_model);
elec_posH=[]; elec_pos_NEW=[];
for i=1:length(elec_comp_h)
    elec_posH(i,:)=elec_comp_h{i}.com; 

    %We can then decompose this to update the direction
    a_i_elec_ii = img_recon.movement_data(i)* ...
        elec_comp_h{i}.tangent(:,1) + ...
        img_recon.movement_data(i+n_elec)* ...
        elec_comp_h{i}.tangent(:,2);  
        
    %Get the old coords of node and update the end point
    elec_pos_NEW(i,1:3) = elec_posH(i,1:3) + a_i_elec_ii'; 
    elec_pos_NEW(i,1:3) = elec_pos_NEW(i,1:3)/norm(elec_pos_NEW(i,1:3));      
end
%Put the normals back in for sphere or clyinder
elec_pos_NEW(:,4:6)=elec_pos_NEW(:,1:3);
mdl_h = ng_mk_gen_models(shape_str, elec_pos_NEW, elec_shape, elec_obj);    

%Create new forward model
img_h=mk_image(mdl_h,1);  
inv_tik_2d.fwd_model=mdl_h;
mdl_h.stimulation= stim; img_h.fwd_model.stimulation=stim;
end

%Make new c2f map and create new images
cfmap = mk_coarse_fine_mapping( mdl_h, mdl_c);
cfmap2 = cfmap./(sum(cfmap,2) * ones(1,size(cfmap,2))); 
img_h=mk_image(mdl_h,1);img_c=mk_image(mdl_c,1); 

%% STEP 2 : SIMULTANEOUS SOLVER

%Get the prior for the conductivity and positions
prior_c=img_c.elem_data; prior_e=zeros(2*n_elec,1);
prior_ce=[prior_c;prior_e];

%Elements on coarse model    
n_elemsc=length(mdl_c.elems(:,1));    
n_elems=length(mdl_h.elems(:,1));
n_elec=length(mdl_h.electrode);
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
mdl_h.jacobian=@jacobian_adjoint;
img_h.fwd_model.jacobian=@jacobian_adjoint;
mdl_h.n_elemsc=n_elemsc; %Assign coarse elements
inv_tik_2d.fwd_model= mdl_h;
inv_tik_2d.fwd_model.coarse2fine = cfmap;
inv_tik_2d.reconst_type='absolute';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.RtR_prior=@prior_movement_tangential;
inv_tik_2d.prior_movement.RegC.func=@prior_laplace;                        
inv_tik_2d.prior_movement.RegM.func=@tikhonov_movement_image_prior;
inv_tik_2d.solve=@inv_solve_abs_GN_diff_tangential;
inv_tik_2d.prior_movement.parameters=hpmt_params;     
inv_tik_2d.hyperparameter.value=hp_params;   
inv_tik_2d.parameters.max_iterations=max_its;

%Add some new fields to the inversion
inv_tik_2d.prior_c   = prior_c;
inv_tik_2d.prior_e   = prior_e;
inv_tik_2d.c2f2 = cfmap2; %Another coarse2 fine for the 
img_h.fwd_model.n_elemsc=n_elemsc;
inv_tik_2d.img_fine=img_h;
inv_tik_2d.mdl_coarse=mdl_c;

%Differene solve on the inverse model (this has mdl_fine and c2f)
img_recon=inv_solve(inv_tik_2d,v_i);
img_reconc=img_c; img_reconc.elem_data=img_recon.elem_data(1:n_elemsc);