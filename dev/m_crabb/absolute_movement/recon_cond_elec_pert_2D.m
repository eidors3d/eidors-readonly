clc; close all; run ~/EIT/Code/mk_paths.m

%Choose absolute solver (1 is Tikhonov and 2 is TV)
solve_type=2;

%Parameter for inclusion conductivity
cond_inc=2; cond_inc2=0.5; move_fac=0.2;

%Pick hyperparameters for the simultaneous solver
hp=10^-4; hpmtohp=8;

%Absolute regularisation parameter
hpabsTik=10^-3; hpabsTV=10^-5; 

%Iterations for absolute movement
move_its=5; cond_its=10;

%No .of electrodes and the anomaly parameter
nelec=16; elec_t=[0.05,0,0.05]; %CEM

%Homogeneous model    TH(ii) = 360/(2*pi)*atan(-x/-y) + 180;
%Create circle homogeneous CEM (0.05) elecs equal spaced
elec_pos=zeros(nelec,2); elec_pos_cart=zeros(nelec,2);
for i=1:nelec
    elec_pos(i,1)=(i-1)*360/nelec; elec_pos(i,2)=0;
end
%Create the homogeneous model    
mdl_h=ng_mk_cyl_models([0,1,0.1],elec_pos,elec_t);    
nodes=mdl_h.nodes;

%Inhomogeneous model
%Create circle inhomogeneous CEM with perturbed electrodes
elec_pos_i=elec_pos;
for i=1:nelec
    elec_pos_i(i,1)= elec_pos(i,1) + move_fac*(-1+2*rand(1,1))*360/nelec;
    elec_pos_i(1,2)=0;
end
%Inhomogeneous model
ball1='solid ball1 = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;';
ball2='solid ball2 = cylinder(-0.2,-0.2,0;-0.2,-0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;';
extra={'ball1','ball2',[ball1,ball2]};
%extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;'};
[mdl_i,mat_idx_i]=ng_mk_cyl_models(0,elec_pos_i,elec_t,extra);

%Make some stimulation patterns and add to models
stim=mk_stim_patterns(16,1,[0 1],[0 1]);
mdl_h.stimulation= stim; 
mdl_i.stimulation= stim; 

%Make inhomogeneous image
img_h=mk_image(mdl_h,1);
img_i=mk_image(mdl_i,1); 
img_i.elem_data(mat_idx_i{2})=cond_inc;
img_i.elem_data(mat_idx_i{3})=cond_inc2;
figure; show_fem(img_h,[1,1,0]);
figure; show_fem(img_i,[1,1,0]);

%Inhomgeneous voltages
v_i=fwd_solve(img_i.fwd_model,img_i);
v_h=fwd_solve(img_h.fwd_model,img_h);     v_hd = v_h.meas; 
figure; hold on; plot(v_h.meas,'b'); plot(v_i.meas,'r'); plot(v_h.meas-v_i.meas,'g'); hold off;

%Calc the inhomogeneous model components
elec_comp_i=calc_electrode_components(img_i.fwd_model);
for i=1:length(elec_comp_i)
   elec_posI(i,:)=elec_comp_i{i}.com;
end

%% NEW MOVEMENT SOLVER
%% SOME EXPERIMENTS AND IDEAS
%1. Movement only then v_h=v_i almost i.e. account for conductivity change
%2. Idea is we have outer iteration. For each of these we do a
%movementandconductivity linearised update. We update the electrode
%positions (assumed moving tangentially) by remeshing. We then do a few
%iterations of absolute, and back to beginning. 

%Outer iteration to do movement every so often
for iii=1:move_its
mdl_h.jacobian=@jacobian_movement_sampling_electrode_tangential;
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
inv_tik_2d.reconst_type='difference';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.fwd_model=mdl_h;
inv_tik_2d.hyperparameter.value=hp; %Appears best hyperparameter
inv_tik_2d.RtR_prior=@prior_electrode_movement;
inv_tik_2d.prior_movement.parameters(1)=hpmtohp;
inv_tik_2d.solve=@inv_solve_diff_GN_one_step;
img_r_tik_tik =inv_solve(inv_tik_2d,v_h,v_i);
n_elems=length(mdl_h.elems(:,1)); 

%Copy and show image
img_r_tik_tikc=img_r_tik_tik;
img_r_tik_tikc.elem_data=img_r_tik_tik.elem_data(1:n_elems);
figure; show_fem(img_r_tik_tikc,[1,1,0]);

%Calculate the electrode components
elec_comp_h=calc_electrode_components(img_h.fwd_model);
for i=1:length(elec_comp_h)
    elec_posH(i,:)=elec_comp_h{i}.com; 
    a_i_elec_ii = img_r_tik_tik.elem_data(n_elems+ i)* ...
        elec_comp_h{i}.tangent;
        
    %Get the old coords of node and update the end point
    elec_pos_NEW(i,:) = elec_posH(i,:) + a_i_elec_ii';          
end
    
%We have updated electrode positions in Cartesians convert to polars
for ii=1:length(elec_pos_NEW(:,1))
    elec_pos_i=elec_pos;
    x=elec_pos_NEW(ii,1);
    y=elec_pos_NEW(ii,2);
    if( x>0 ) 
        if(y>0)
            TH(ii) = 360/(2*pi)*atan(x/y);
        elseif(y<0)
            TH(ii) = 360/(2*pi)*atan(-y/x) + 90;
        end
    elseif( x < 0)
        if( y > 0 )
            TH(ii) = 360/(2*pi)*atan(y/-x) + 270;
        elseif(y<0)            
            TH(ii) = 360/(2*pi)*atan(-x/-y) + 180;
        end
    end
    elec_pos(ii,1)=TH(ii);
end
    
%Remesh the model and then create coarse fine mapping with image
mdl_h=ng_mk_cyl_models([0,1,0.1],elec_pos,elec_t);    
mdl_h.stimulation= stim; 
    
%Make coarse to fine make
c2f=mk_coarse_fine_mapping(mdl_h,img_r_tik_tik.fwd_model);
Nelem_data = c2f*(img_h.elem_data + img_r_tik_tik.elem_data(1:n_elems)); 

%Create new image for iterations
img_h=mk_image(mdl_h,Nelem_data);  
figure; show_fem(img_h,[1,1,0]);

%Now peform an absolute solve
mdl_h.jacobian=@jacobian_adjoint;
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
inv_tik_2d.reconst_type='absolute';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.fwd_model=mdl_h;

if(solve_type==1)
    inv_tik_2d.RtR_prior=@prior_laplace;
    inv_tik_2d.hyperparameter.value=hpabsTik; %Appears best hyperparameter
    inv_tik_2d.solve=@inv_solve_abs_GN;        
elseif(solve_type==2)
    inv_tik_2d.R_prior=@prior_TV;
    inv_tik_2d.hyperparameter.value=hpabsTV; %Appears best hyperparameter
    inv_tik_2d.solve=@inv_solve_abs_pdipm;    
end    
inv_tik_2d.parameters.max_iterations=cond_its;
inv_tik_2d.parameters.term_tolerance=10^-10;
img_r_tik_tik =inv_solve(inv_tik_2d,v_i);
n_elems=length(mdl_h.elems(:,1)); 

figure; show_fem(img_r_tik_tik,[1,1,0]);

%Reiterate by creating new model and image
mdl_h=img_r_tik_tik.fwd_model;
img_h=mk_image(mdl_h,img_r_tik_tik.elem_data(1:n_elems));
v_h=fwd_solve(img_h);

end

%Plot the final
figure; hold on; plot(v_h.meas,'b'); plot(v_i.meas,'r'); plot(v_h.meas-v_i.meas,'g'); hold off;
