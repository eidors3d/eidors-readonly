%clear all; close all; run ~/EIT/Code/mk_paths.m; clc;

%Conductivity and permittivity parametersma
cond_back=1; cond_inhom=100; 
perm_back=0.01; perm_inhom=0.01;
%perm_back=0; perm_inhom=0;

%Number of Fourier coefficients and time interval [-T_max,T_max]
T_max = 0.5; n_freq = 1001;
%Plotting interval [T_min_plot,T_max_plot] with T_freq samples
T_min_plot = -0.02; T_max_plot = 0.08; T_freq = 700; %700<25MB gmail image!!

%Contact resistance and capacitance parameters
R_s = 0.01;
R_con = 0.01; 
C_con_list = [0];

%Sweep through frequencies (NB odd coeffs only for square wave)
n_freqs =-n_freq:2:n_freq;
T_freqs = linspace(T_min_plot,T_max_plot,T_freq);

%Make forward model
ball1 = 'solid ball = cylinder(0.3,0.3,0;0.3,0.3,1;0.3) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;';
box1  = 'solid box = orthobrick(-0.4,-0.4,0;-0.1,0.1,0.05) -maxh=0.1;';
extra = {'ball','box',[ball1,box1]};
fmdl= ng_mk_cyl_models(0,[8],[0.2,0,0.05],extra); 
figure; show_fem(fmdl);

%Make stimulation pattern
sp_mp = [1, 5, 2, 3]; %+I,-I,+V,-V
stim = stim_meas_list(sp_mp,8,.01,1);
fmdl.stimulation = stim;

%Storage for the potentials
volt_nodes_omega_freqs = zeros(size(fmdl.nodes,1),length(n_freqs),length(C_con_list));
volt_elecs_omega_freqs = zeros(1,length(n_freqs),length(C_con_list));
volt_nodes_T_freqs     = zeros(size(fmdl.nodes,1),length(T_freqs),length(C_con_list));
volt_elecs_T_freqs     = zeros(1,length(T_freqs),length(C_con_list));
I_in=zeros(1,length(T_freqs),1);
f_k=zeros(2*n_freq+1);

%Loop through some capacitances
for ccc=1:length(C_con_list) 
%Get capacitance
C_con=C_con_list(ccc);
fprintf(1,sprintf('Capacitance no. %6.2f of %6.2f\n',ccc,length(C_con_list)));

%Solve PDE for integers n: \nabla \cdot( (\sigma + \epsilon i n \pi / T_max) \nabla u) 
for kk=1:length(n_freqs)
%Angular frequency
kk_freq = n_freqs(kk); %angular frequency
omega_freq = kk_freq*pi/T_max;
fprintf(1,sprintf('Frequency no. %6.2f of %6.2f\n',kk,length(n_freqs)));

%Add complex contact impedance as function of frequency
for ff=1:size(fmdl.electrode,2)
   fmdl.electrode(ff).z_contact = R_s + 1/( 1/R_con + 1i*omega_freq*C_con ); 
   fmdl.electrode(ff).z_contact = 1/( 1/R_con + 1i*omega_freq*C_con );    
end

%Make image which is frequency dependent (assume linear)
img = mk_image( fmdl, cond_back + 1i*perm_back*omega_freq );
img.elem_data( fmdl.mat_idx{2} ) = cond_inhom + 1i*perm_back*omega_freq;
img.elem_data( fmdl.mat_idx{3} ) = cond_back  + 1i*perm_inhom*omega_freq;
img.fwd_solve.get_all_meas=1; %Get all measurements

%Plot the target real and complex
if(0)
figure;
img.calc_colours.component = 'real';
subplot(121); show_fem(img,[1,0,0]);
img.calc_colours.component = 'imag';
subplot(122); show_fem(img,[1,0,0]);
end

%Compute Fourier coefficient and put as current RHS
f_k(kk+n_freq+1) = -2/(1i*kk_freq*pi);
stim.stim_pattern(1,1) = 0.01*f_k(kk+n_freq+1);
stim.stim_pattern(5,1) = -0.01*f_k(kk+n_freq+1);
fmdl.stimulation=stim; img.fwd_model.stimulation=stim;

%Forward solve for voltage and and add to lists
vi = fwd_solve(img);
volt_nodes_omega_freqs(:,kk,ccc) = vi.volt(:,1);
volt_elecs_omega_freqs(:,kk,ccc) = vi.meas; %Single electrode pait

end

%Compute the interior/electrode potential and current pattern
for jj=1:length(T_freqs)
    fprintf(1,sprintf('Time no. %6.2f of %6.2f\n',jj,length(T_freqs)));
    jj_T = T_freqs(jj);
    for kk = 1:length(n_freqs)
        kk_freq=n_freqs(kk);
        omega_freq = kk_freq*pi/T_max;        
        volt_nodes_T_freqs(:,jj,ccc) = volt_nodes_T_freqs(:,jj,ccc) + volt_nodes_omega_freqs(:,kk,ccc)*exp(1i*omega_freq*jj_T);
        volt_elecs_T_freqs(:,jj,ccc) = volt_elecs_T_freqs(:,jj,ccc) + volt_elecs_omega_freqs(:,kk,ccc)*exp(1i*omega_freq*jj_T);    
        I_in(jj) = I_in(jj) + f_k(kk+n_freq+1)*exp(1i*omega_freq*jj_T) ;         
    end
  
end

%Movie variable colorbar
figure; imgn = rmfield(img,'elem_data');
for jj=1:length(T_freqs)
    fprintf(1,sprintf('Time no. %6.2f of %6.2f\n',jj,length(T_freqs)));
    jj_T = T_freqs(jj);    
    imgn.node_data = volt_nodes_T_freqs(:,jj,ccc);
    %figure;
    imgn.calc_colours.component = 'real';
    show_fem(imgn,[1,0,0]);
    title(sprintf('Potential at T= %1.3f ms and C_con = %7.6f and Perhom_inhom = %7.6f',jj_T,C_con,perm_inhom));   
    M(jj) = getframe(gcf);      
end

%Plot electrode current as function of time
figure; hold on;
plot(T_freqs,I_in,'r');
legend('Current: E1-E5');
title(sprintf('Electrode input current as function of time C_con = %7.6f and Perhom_inhom = %7.6f',C_con,perm_inhom));
saveas(gcf,sprintf('electrode_currents_square_C_con_%7.6f_perhom_inhom_%7.6f.eps',C_con,perm_inhom),'epsc');

%Plot electrode potential as function of time
figure; hold on;
plot(T_freqs,real(volt_elecs_T_freqs(1,:,ccc)),'r*');
legend('Measure: E2-E3');
title(sprintf('Electrode measured potential as function of time C_con = %7.6f and Perhom_inhom = %7.6f',C_con,perm_inhom));
saveas(gcf,sprintf('electrode_voltages_square_C_con_%7.6f_perhom_inhom_%7.6f.eps',C_con,perm_inhom),'epsc');

end
