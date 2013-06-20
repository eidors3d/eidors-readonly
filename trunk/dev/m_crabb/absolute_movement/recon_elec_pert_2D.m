clc; close all; run ~/EIT/Code/mk_paths.m

%% PARAMETERS
%Parameter for inclusion conductivity
cond_inc=1; move_fac=0.3; elec_ref=0.01; dom_ref=0.05;
%No .of electrodes and the anomaly parameter
nelec=16; elec_t=[0.05,0,elec_ref]; %CEM

%% Homogeneous model    TH(ii) = 360/(2*pi)*atan(-x/-y) + 180;
%Create circle homogeneous CEM (0.05) elecs equal spaced
elec_pos=zeros(nelec,2); elec_pos_cart=zeros(nelec,2);
for i=1:nelec
    elec_pos(i,1)=(i-1)*360/nelec; %Angle
    elec_pos(i,2)=0; %Height
end
%Create the homogeneous model    
mdl_h=ng_mk_cyl_models([0,1,dom_ref],elec_pos,elec_t);    
nodes=mdl_h.nodes;

%% Inhomogeneous model
%Create circle inhomogeneous CEM with perturbed electrodes
elec_pos_i=elec_pos;
for i=1:nelec
    elec_pos_i(i,1)= elec_pos(i,1) + move_fac*(-1+2*rand(1,1))*360/nelec; %Angle
    elec_pos_i(1,2)=0; %Height
end
%Inhomogeneous model
extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;'};
[mdl_i,mat_idx_i]=ng_mk_cyl_models([0,1,dom_ref],elec_pos_i,elec_t,extra);

%Make some stimulation patterns and add to models
stim=mk_stim_patterns(16,1,[0 1],[0 1]);
mdl_h.stimulation= stim; mdl_i.stimulation= stim; 

%Make inhomogeneous image
img_h=mk_image(mdl_h,1);
img_i=mk_image(mdl_i,1); img_i.elem_data(mat_idx_i{2})=cond_inc;
figure; show_fem(img_h,[1,1,0]);
figure; show_fem(img_i,[1,1,0]);

%Inhomgeneous voltages
v_i=fwd_solve(img_i.fwd_model,img_i);
v_h=fwd_solve(img_h.fwd_model,img_h);  v_hd = v_h.meas; 

%Calc the inhomogeneous model components
elec_comp_i=calc_electrode_components(img_i.fwd_model);
elec_comp=calc_electrode_components(img_h.fwd_model);

%True positions (to compute some errors)
for i=1:length(elec_comp_i)
   elec_posI(i,:)=elec_comp_i{i}.com;
end
%Starting positions (to iteratively update)
for i=1:length(elec_comp)
   elec_posH(i,:)=elec_comp{i}.com; 
end
elec_posINITIAL=elec_posH;

%Initial error
error_pos(1) = norm(elec_posH-elec_posI);    
error_volt(1) = norm(v_hd-v_i.meas);    

%Regularisation parameter and iterations
alpha=0.00001; its=4; 

%Loop through the iterations and perform Armijo GN
for k=1:its        

%Movement Jacobian, (Tikhonov) regularisation matrix and parameter
J=jacobian_electrode_movement_sampling(img_h);
reg_dim=size(J'*J,1); reg=eye(reg_dim); 

%Update from GN (do Armijo linesearch??)
GN_update = (J'*J + alpha*reg)\(-J'*(v_hd-v_i.meas));

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
mdl_h=ng_mk_cyl_models([0,1,dom_ref],elec_pos,elec_t); nodes=mdl_h.nodes;
mdl_h.stimulation= stim; 
img_h=mk_image(mdl_h,1);
figure; show_fem(img_h,[1,1,0]);

%Get the simulated voltages
v_h=fwd_solve(img_h.fwd_model,img_h);    
v_hd = v_h.meas; 

%Calculate update electrode configuration
elec_comp=calc_electrode_components(img_h.fwd_model);
for i=1:length(elec_comp)
    elec_posH(i,:)=elec_comp{i}.com;
end

%The error in the positions from actual
error_pos(k+1,1) = norm(elec_posH-elec_posI);    
error_volt(k+1,1) = norm(v_hd-v_i.meas);    
    
end

%Plot the erros in coordinates
figure; plot(error_volt);
