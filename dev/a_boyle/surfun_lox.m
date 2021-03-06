function [sigma_true,sigma_ref,sigma_rec,errm,irrm,mshd,msht] = surfun_lox(tfac,p,t,k)
%Function that formulates and solves the inverse problem in ERT using the 
%LOG BOUNDS surrogate functions, on a 2d square grid with point electrodes
%
%INPUTS
% tfac - Tikhonov factor, p - lower bound for k * log conductivity, t - upper
% bound for k * log conductivity, k - kappa the surrogate scaling factor
%
% Target conductivity below is between 0.1 and 10.921 
%
%Example call
%[sigma_true,sigma_ref,sigma_rec,errm,irrm,mshd,msht] =
%surfun_lox(1e-6,-3,3,1);
%
%To plot the reconstructions use
%figure; subplot(1,3,1); plotinvsol(sigma_true,mshd.vtx,mshd.simp);
%subplot(1,3,2); plotinvsol(sigma_rec(:,1),msht.vtx,msht.simp);
%subplot(1,3,3); plotinvsol(sigma_rec(:,2),msht.vtx,msht.simp);


%Load forward mesh, fine grid with 30 point electrodes presaved
load regsquare_point mshd %(data mesh - fine resolution)
%load regsquare_point msht %(solution mesh - coarse resolution)

vtx=mshd.vtx; simp=mshd.simp; elec=mshd.elec; gnd_ind=mshd.gnd_ind;

%The centers of the elements
cc = zeros(size(simp,1),size(vtx,2));
for ii=1:length(simp)
    cc(ii,:) = mean(vtx(simp(ii,:),:));
end


%Choose a target conductivity distribution 
sigma_true = 1 - 0.9*exp(-((cc(:,1)-(-5)).^2 + (cc(:,2) - (-24)).^2)/(2*(5).^2))+...
    10*exp(-((cc(:,1)-5).^2 + (cc(:,2) - (-10)).^2)/(2*3.^2))+...
    8*exp(-((cc(:,1)-(-8)).^2 + (cc(:,2) - (-8)).^2)/(2*3.^2))+...
    6*exp(-((cc(:,1)-(10)).^2 + (cc(:,2) - (-20)).^2)/(2*3.^2));


%Define rho_true based on sigma_true
rho_true = -k*log(((t-p)/k)*ones(size(simp,1),1)./(log(sigma_true) - ...
    (p/k)*ones(size(simp,1),1)) - ones(size(simp,1),1));


disp([exp(p/k), min(sigma_true), max(sigma_true), exp(t/k)])

%Set boundary currents
protocol = '{op}'; scheme = '{normal}';
[I] = set_currents(protocol,elec,vtx,gnd_ind);
[E,~] = fem_master(vtx,simp,sigma_true,scheme,gnd_ind);
V=E\I; [voltage] = voltage_meas(protocol,V,elec);

%Add noise to simulate the data 
Ce = 1e-2*abs(mean(voltage))*eye(length(voltage));

%The "experimental data" 
data = voltage + sqrt(Ce)*randn(length(voltage),1);

%Import model for coarse mesh
load regsquare_point msht %(solution mesh - coarse resolution)
vtx=msht.vtx; simp=msht.simp; elec=msht.elec; gnd_ind=msht.gnd_ind;

%True coarse images
cc = zeros(size(simp,1),size(vtx,2));
for ii=1:length(simp)
    cc(ii,:) = mean(vtx(simp(ii,:),:));
end

sigma_true_co = 1 -0.9*exp(-((cc(:,1)-(-5)).^2 + (cc(:,2) - (-24)).^2)/(2*(5).^2))+...
    10*exp(-((cc(:,1)-5).^2 + (cc(:,2) - (-10)).^2)/(2*3.^2))+...
    8*exp(-((cc(:,1)-(-8)).^2 + (cc(:,2) - (-8)).^2)/(2*3.^2))+...
    6*exp(-((cc(:,1)-(10)).^2 + (cc(:,2) - (-20)).^2)/(2*3.^2));


rho_true_co = -k*log(((t-p)/k)*ones(size(simp,1),1)./(log(sigma_true_co) - ...
    (p/k)*ones(size(simp,1),1)) - ones(size(simp,1),1));


%Check mats for consistency 
if abs(mean(sigma_true)-mean(sigma_true_co)) > 1 
    warning('Missmatched targets')
end

%Regularization details
%L = iso_f_smooth(simp,vtx,1,1); 
L = eye(size(simp,1));


%Initialize at homogeneous estimate
[I] = set_currents(protocol,elec,vtx,gnd_ind);

sigma_ref = ((exp(p/k) + exp(t/k))/2)*ones(length(simp),1);

rho_ref = -k*log(((t-p)/k)*ones(size(simp,1),1)./(log(sigma_ref) - ...
    (p/k)*ones(size(simp,1),1)) - ones(size(simp,1),1));


%Initialize sol ans stat vectors
sigma_rec = zeros(size(simp,1),2); 
errm = zeros(2,1); 
irrm = zeros(2,1); 


%% First Iteration 

%Solve the forward problem and compute
[Eref,~] = fem_master(vtx,simp,sigma_ref,scheme,gnd_ind);
Vref = Eref\I; [ref] = voltage_meas(protocol,Vref,elec);


%Voltage residual
dva(:,1) = data - ref;


%The Jacobians
J = J_pe(vtx,simp,elec,length(data(:,1)),sigma_ref,gnd_ind,protocol);
nu_m = exp(p/k) * exp((ones(size(simp,1),1)*(t-p)/k)./...
    (ones(size(simp,1),1) + exp(-rho_ref/k)));
nu_b_minus_p = (t-p)*ones(size(simp,1),1)./(ones(size(simp,1),1) + exp(-rho_ref/k));
nu1_minus_m = ones(size(simp,1),1)./(ones(size(simp,1),1) + exp(rho_ref/k));
Delta = diag((nu_m.*nu1_minus_m.*nu_b_minus_p)*(k^-2));
A = J*Delta;


%Solve surrogate problem for m with tfacm
%Record the solution sigma_rec(:,it)
rho_rec(:,1) = rho_ref + (A.'*A + tfac*L.'*L)\(A.'*dva(:,1));
sigma_rec(:,1) = (exp(p/k)* ones(size(simp,1),1).*exp(((t-p)/k)./(ones(size(simp,1),1) ...
+ exp(-rho_rec(:,1)/k))) );


%Record its data error err(:,it)
[Eerrm,~] = fem_master(vtx,simp,sigma_rec(:,1),scheme,gnd_ind);
Verrm = Eerrm\I; [ref_errm] = voltage_meas(protocol,Verrm,elec);
errm(1) = norm(data - ref_errm)/norm(data);

%Record its image error irr(:,it)
irrm(1) = norm(sigma_true_co - sigma_rec(:,1))/norm(sigma_true_co);

%% Second iteration

%Solve the surrogate problem's second iteration
[Eref,~] = fem_master(vtx,simp,sigma_rec(:,1),scheme,gnd_ind);
Vref = Eref\I; [ref] = voltage_meas(protocol,Vref,elec);

%Voltage residual
dva(:,2) = data - ref;

%The Jacobian for sigma
J = J_pe(vtx,simp,elec,length(data),sigma_rec(:,1),gnd_ind,protocol);
nu_m = exp(p/k) * exp((ones(size(simp,1),1)*(t-p)/k)./...
    (ones(size(simp,1),1) + exp(-rho_rec(:,1)/k)));
nu_b_minus_p = (t-p)*ones(size(simp,1),1)./(ones(size(simp,1),1) + exp(-rho_rec(:,1)/k));
nu1_minus_m = ones(size(simp,1),1)./(ones(size(simp,1),1) + exp(rho_rec(:,1)/k));
Delta = diag((nu_m.*nu1_minus_m.*nu_b_minus_p)*(k^-2));
A = J*Delta;

rho_rec(:,2) = rho_rec(:,1) + (A.'*A + tfac*L.'*L)\(A.'*dva(:,2));
sigma_rec(:,2) = (exp(p/k)* ones(size(simp,1),1).*exp(((t-p)/k)./(ones(size(simp,1),1) ...
+ exp(-rho_rec(:,2)/k))) );


%Record its data error err(:,it)
[Eerrm,~] = fem_master(vtx,simp,sigma_rec(:,2),scheme,gnd_ind);
Verrm = Eerrm\I; [ref_errm] = voltage_meas(protocol,Verrm,elec);
errm(2) = norm(data - ref_errm)/norm(data);

%Record its image error irr(:,it)
irrm(2) = norm(sigma_true_co - sigma_rec(:,2))/norm(sigma_true_co);