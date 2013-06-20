clc; close all; run ~/EIT/Code/mk_paths.m

%% PARAMETERS
%Parameter for inclusion conductivity
cond_inc=2; move_fac=0.3; elec_ref=0.05; dom_ref=0.5; dom_h=1;
%No .of electrodes and the anomaly parameter
nelec=16; elec_t=[0.05,0,elec_ref]; %CEM

%Pick the difference hyperparm
hp=10^-4; hpmtohp=5; hpabs=10^-4; 
%Iterations for absolute movement
move_its=10; cond_its=6;

%Homogeneous model    TH(ii) = 360/(2*pi)*atan(-x/-y) + 180;
%Create circle homogeneous CEM (0.05) elecs equal spaced
elec_pos=zeros(nelec,2); elec_pos_cart=zeros(nelec,2);
for i=1:nelec
    elec_pos(i,1)=(i-1)*360/nelec; 
    elec_pos(i,2)=0.5;
end
%Create the homogeneous model    
mdl_h=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos,elec_t);    
nodes=mdl_h.nodes;

%Inhomogeneous model
%Create circle inhomogeneous CEM with perturbed electrodes
elec_pos_i=elec_pos;
for i=1:nelec
    elec_pos_i(i,1)= elec_pos(i,1) + move_fac*(-1+2*rand(1,1))*360/nelec;
    elec_pos_i(i,2) = 0.5 + 0.1*(-1+2*rand(1,1)); %Height
end
%Inhomogeneous model
extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0.2;1,1,0.8) -maxh=0.1;'};
[mdl_i,mat_idx_i]=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos_i,elec_t,extra);

%Make some stimulation patterns and add to models
stim=mk_stim_patterns(16,1,[0 1],[0 1]);
mdl_h.stimulation= stim; 
mdl_i.stimulation= stim; 

%Make inhomogeneous image
img_h=mk_image(mdl_h,1);
img_i=mk_image(mdl_i,1); img_i.elem_data(mat_idx_i{2})=cond_inc;
figure; show_fem(img_h,[1,1,0]);
figure; show_fem(img_i,[1,1,0]);

%Inhomgeneous voltages
v_i=fwd_solve(img_i.fwd_model,img_i);
v_h=fwd_solve(img_h.fwd_model,img_h);     v_hd = v_h.meas; 

%Calc the inhomogeneous model components
elec_comp_i=calc_electrode_components(img_i.fwd_model);
for i=1:length(elec_comp_i)
   elec_posI(i,:)=elec_comp_i{i}.com;
end

%% LINEARISED UPDATE WITH NEW POSITIONS
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
inv_tik_2d.reconst_type='difference';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.fwd_model=mdl_h;
inv_tik_2d.hyperparameter.value=hp; %Appears best hyperparameter
inv_tik_2d.RtR_prior=@prior_tikhonov;
inv_tik_2d.solve=@inv_solve_diff_GN_one_step;
img_r_tik_tik =inv_solve(inv_tik_2d,v_h,v_i);
figure; show_fem(img_r_tik_tik,[1,1,0]);

%% EIDORS MOVEMENT SOLVER
mdl_h.jacobian=@jacobian_movement;
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
inv_tik_2d.reconst_type='difference';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.fwd_model=mdl_h;
inv_tik_2d.hyperparameter.value=hp; %Appears best hyperparameter
inv_tik_2d.RtR_prior=@prior_movement;
inv_tik_2d.prior_movement.parameters(1)=hpmtohp;
inv_tik_2d.solve=@inv_solve_diff_GN_one_step;
img_r_tik_tik =inv_solve(inv_tik_2d,v_h,v_i);
n_elems=length(mdl_h.elems(:,1));    

figure; show_fem_moveCOPY(img_r_tik_tik);
img_r_tik_tik.elem_data=img_r_tik_tik.elem_data(1:n_elems);
figure; show_fem(img_r_tik_tik,[1,1,0]);


%% NEW MOVEMENT SOLVER
%% SOME EXPERIMENTS AND IDEAS
%1. IF WE HAVE FIRST STAGE MOVEMENT ONLY, THEN WE GET V_H=V_I ALMOST AND
%THEN ITERATIVE SOLVER THINKS CONVERGENCE. 
%2. IF WE HAVE SIMULTANEOUS CONDUCTIVITY AND MOVEMENT THEN THIS DOES NOT
%HAPPEN. WE WANT INITIAL ITERATIONS (AS PARAMETER) SO WE PEFORM SOME
%SIMULTANEOUS CONDUCTIVITY MOVEMENT SOLUTION
%3. AFTER (SUFFICIENTLY) MANY INITIAL ITERATIONS WE WANT TO BREAK INTO AN
%ABSOLUTE SOLVER, WHEN THE MOVEMENT IS NOT AN ISSUE ANYMORE. 
%4. CAN I LET GO OF MOVEMENT THOUGH OR DOES THIS HAVE TO KEEP ON BEING
%ACCOUNTED FOR??????

%Outer iteration to do movement every so often
for iii=1:move_its
mdl_h.jacobian=@jacobian_movement_eidors_electrode_tangential;
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

    %We can then decompose this to update the direction
    a_i_elec_ii = img_r_tik_tik.elem_data(n_elems+i)*elec_comp_h{i}.tangent(:,1) + ...
                  img_r_tik_tik.elem_data(n_elems+i+nelec)*elec_comp_h{i}.tangent(:,2);    
        
    %Get the old coords of node and update the end point
    elec_pos_NEW(i,:) = elec_posH(i,:) + a_i_elec_ii';          
end
    
%We have updated electrode positions in Cartesians convert to polars
for ii=1:length(elec_pos_NEW(:,1))
    elec_pos_i=elec_pos;
    x=elec_pos_NEW(ii,1);
    y=elec_pos_NEW(ii,2);
    z=elec_pos_NEW(ii,3);
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
    elec_pos(ii,2)=z;
end
    
%Remesh the model and then create coarse fine mapping with image
mdl_h=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos,elec_t);    
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
inv_tik_2d.hyperparameter.value=hpabs; %Appears best hyperparameter
inv_tik_2d.RtR_prior=@prior_laplace;
%inv_tik_2d.intial_c=img_h.elem_data; %CHANGE BFGS SOLVER!!!!
%inv_tik_2d.prior_c=img_h.elem_data;
inv_tik_2d.solve=@inv_solve_abs_BFGS;
inv_tik_2d.parameters.max_iterations=cond_its;
inv_tik_2d.parameters.term_tolerance=10^-12;
img_r_tik_tik =inv_solve(inv_tik_2d,v_i);
n_elems=length(mdl_h.elems(:,1)); 

figure; show_fem(img_r_tik_tik,[1,1,0]);

%Reiterate by creating new model and image
mdl_h=img_r_tik_tik.fwd_model;
img_h=mk_image(mdl_h,img_r_tik_tik.elem_data(1:n_elems));
v_h=fwd_solve(img_h);

end

%% LINEARISED UPDATE WITH NEW POSITIONS
%{
its=1;
for iii=1:its
n_elems=length(mdl_h.elems(:,1));    
%Calc the homogeneous components
elec_comp_h=calc_electrode_components(img_h.fwd_model);
for i=1:length(elec_comp_h)
   elec_posH(i,:)=elec_comp_h{i}.com; 
end
v_h=fwd_solve(img_h.fwd_model,img_h);  v_hd = v_h.meas; 

%Initialise an inverse model
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
inv_tik_2d.reconst_type='difference';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.fwd_model=mdl_h;
inv_tik_2d.hyperparameter.value=10^-4;
inv_tik_2d.prior_movement.parameters(1)=5;
inv_tik_2d.RtR_prior=@prior_laplace;

img_r_tik_tik_elec = inv_solve_diff_GN_one_step_movement(inv_tik_2d,v_h,v_i);

%Loop through the electrode components
cnt=0;
for i=1:length(elec_comp_h); %for each electrode
    %We can then decompose this to update the direction
    a_i_elec_ii = img_r_tik_tik_elec.elem_data(n_elems+ i)* ...
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
c2f=mk_coarse_fine_mapping(mdl_h,img_r_tik_tik_elec.fwd_model);
Nelem_data = c2f*(img_h.elem_data + img_r_tik_tik_elec.elem_data(1:n_elems)); 

%Create new image for iterations
img_h=mk_image(mdl_h,Nelem_data);  
figure; show_fem(img_h,[1,1,0]);
            
end
%}
%% ITERATIVELY UPDATE ELECTRODE POSITIONS
%INITIALLY UPDATING ELECTRODES WITHOUT CONDUCTIVITY IS BAD!!!
%VOLTAGE MISMATCH IS EFFECTIVELY ZERO FOR THIS
%{
if(solve_type==2)
%Regularisation parameter and iterations
alpha=0.0001; its=5;
   
%STEP 1 - INITIAL MOVEMENT CORRECTION   
%THIS GIVES ZERO ERROR IN THE VOLTAGES!!!!!!!!!!
%Movement Jacobian, (Tikhonov) regularisation matrix and parameter
J=jacobian_electrode_movement_sampling(img_h);
reg_dim=size(J'*J,1); reg=eye(reg_dim); 

%Update from GN (do Armijo linesearch??)
GN_update = (J'*J + alpha*reg)\(-J'*(v_hd-v_i.meas));

%Calc the inhomogeneous model components
elec_comp=calc_electrode_components(img_h.fwd_model);
for i=1:length(elec_comp)
   elec_posH(i,:)=elec_comp{i}.com;
end

%Update the COM position from avergae tangential movement
for i=1:length(elec_comp);           
    %We can then decompose this to update the direction
    a_i_elec_ii = GN_update(i)*elec_comp{i}.tangent;
        
    %Get the old coords of node and update the end point
    elec_pos_NEW(i,:) = elec_posH(i,:) + a_i_elec_ii';     
end

%Updated electrode positions in Cartesians so convert to Polars
for ii=1:length(elec_pos_NEW(:,1))
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

%Create the homogeneous model    
mdl_h=ng_mk_cyl_models([0,1,0.1],elec_pos,elec_t); nodes=mdl_h.nodes;
mdl_h.stimulation= stim; 
img_h=mk_image(mdl_h,1);
figure; show_fem(img_h,[1,1,0]);

%Get the simulated voltages
v_h=fwd_solve(img_h.fwd_model,img_h);    
v_hd = v_h.meas; %Initialise an inverse model

%STEP 2 - ABSOLUTE CONDUCTIVITY
inv_tik_2d=eidors_obj('inv_model','EIT inverse');
inv_tik_2d.reconst_type='absolute';
inv_tik_2d.jacobian_bkgnd.value=img_h.elem_data;
inv_tik_2d.fwd_model=mdl_h;
inv_tik_2d.hyperparameter.value=5*10^-4;
inv_tik_2d.solve=@inv_solve_abs_GN;
inv_tik_2d.RtR_prior=@prior_tikhonov;
inv_tik_2d.parameter.max_iterations=its;

%Solve
img=inv_solve(inv_tik_2d,v_i);
figure; show_fem(img,[1,1,0]);

end
%}
