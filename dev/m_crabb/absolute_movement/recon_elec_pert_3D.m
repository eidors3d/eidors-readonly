clc; close all; run ~/EIT/Code/mk_paths.m

%% PARAMETERS
%Parameter for inclusion conductivity
cond_inc=1; move_fac=0.3; elec_ref=0.05; dom_ref=0.1; dom_h=1;
%No .of electrodes and the anomaly parameter
nelec=16; elec_t=[0.05,0,elec_ref]; %CEM
%Regularisation parameter and iterations
alpha=0.0001; its=8; 

%% Homogeneous model    TH(ii) = 360/(2*pi)*atan(-x/-y) + 180;
%Create circle homogeneous CEM (0.05) elecs equal spaced
elec_pos=zeros(nelec,2); elec_pos_cart=zeros(nelec,2);
for i=1:nelec
    elec_pos(i,1)=(i-1)*360/nelec; %Angle
    elec_pos(i,2)=0.5; %Height
end
%Create the homogeneous model    
mdl_h=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos,elec_t);    
nodes=mdl_h.nodes;

%% Inhomogeneous model
%Create circle inhomogeneous CEM with perturbed electrodes
elec_pos_i=elec_pos;
for i=1:nelec
    elec_pos_i(i,1) = elec_pos(i,1) + move_fac*(-1+2*rand(1,1))*360/nelec; %Angle
    elec_pos_i(i,2) = 0.5 + 0.1*(-1+2*rand(1,1)); %Height
end
%Inhomogeneous model
extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;'};
[mdl_i,mat_idx_i]=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos_i,elec_t,extra);

%% Inhomogeneous model (pertrub theta but not z)
elec_pos_i_non_z=elec_pos_i;
for i=1:nelec
    elec_pos_i_non_z(i,2)=0.5;    
end
%Inhomogeneous model
extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;'};
[mdl_i_non_z,mat_idx_i_non_z]=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos_i_non_z,elec_t,extra);


%Make some stimulation patterns and add to models
stim=mk_stim_patterns(16,1,[0 1],[0 1]);
mdl_h.stimulation= stim; mdl_i.stimulation= stim; mdl_i_non_z.stimulation=stim;

%Make inhomogeneous image
img_h=mk_image(mdl_h,1);
img_i=mk_image(mdl_i,1); img_i.elem_data(mat_idx_i{2})=cond_inc;
img_i_non_z=mk_image(mdl_i_non_z,1); img_i_non_z.elem_data(mat_idx_i_non_z{2})=cond_inc;
figure; show_fem(img_h,[1,1,0]);
figure; show_fem(img_i,[1,1,0]);
figure; show_fem(img_i_non_z,[1,1,0]);

%Inhomgeneous voltages
v_i=fwd_solve(img_i.fwd_model,img_i);
v_h=fwd_solve(img_h.fwd_model,img_h);  v_hd = v_h.meas; 
v_i_non_z=fwd_solve(img_i_non_z.fwd_model,img_i_non_z);
figure; plot(v_h.meas,'r'), hold on; plot(v_i.meas,'g'); hold on; plot(v_i_non_z.meas,'b'); hold off;

%Calc the inhomogeneous model components
elec_comp_i=calc_electrode_components(img_i.fwd_model);
elec_comp=calc_electrode_components(img_h.fwd_model);

%Plot the normals and tangents of each electrode
%{
figure; show_fem(mdl_h,[1,1,0]); hold on;
for i=1:length(elec_comp)
    com=elec_comp{i}.com;
    normal=elec_comp{i}.normal;
    tangential=elec_comp{i}.tangent;
    quiver3(com(1),com(2),com(3),normal(1),normal(2),normal(3),'r*');
    quiver3(com(1),com(2),com(3),tangential(1,1),tangential(2,1),tangential(3,1),'b*');
    quiver3(com(1),com(2),com(3),tangential(1,2),tangential(2,2),tangential(3,2),'b*');    
end
hold off
%}

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

%Loop through the iterations and perform Armijo GN
for k=1:its        

%Movement Jacobian, (Tikhonov) regularisation matrix and parameter
%J=jacobian_electrode_movement_sampling(img_h);
J=jacobian_movement_eidors_electrode_tangential(img_h);
n_elem=length(img_h.fwd_model.elems(:,1)); J = J(:,n_elem+1:end);
reg_dim=size(J'*J,1); reg=eye(reg_dim); 

%Update from GN (do Armijo linesearch??)
GN_update = (J'*J + alpha*reg)\(-J'*(v_hd-v_i.meas));

%Update the COM position from avergae tangential movement
for i=1:length(elec_comp);           
    %We can then decompose this to update the direction
    a_i_elec_ii = GN_update(i)*elec_comp{i}.tangent(:,1) + ...
                  GN_update(i+nelec)*elec_comp{i}.tangent(:,2);
        
    %Get the old coords of node and update the end point
    elec_pos_NEW(i,:) = elec_posH(i,:) + a_i_elec_ii';     
end

%Updated electrode positions in Cartesians so convert to Polars
for ii=1:length(elec_pos_NEW(:,1))
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
    %Put into the cylindrical coordinates
    elec_pos(ii,1)=TH(ii);
    elec_pos(ii,2)=z;
end

%Create the homogeneous model    
mdl_h=ng_mk_cyl_models([dom_h,1,dom_ref],elec_pos,elec_t); nodes=mdl_h.nodes;
mdl_h.stimulation= stim; 
img_h=mk_image(mdl_h,1);
figure; show_fem(img_h,[1,1,0]);

%Get the simulated voltages
v_h=fwd_solve(img_h.fwd_model,img_h);    
v_hd = v_h.meas; 

%Plot the simulated and measured voltages
figure; plot(v_h.meas,'r'), hold on; plot(v_i.meas,'g');

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
figure; plot(error_pos);
