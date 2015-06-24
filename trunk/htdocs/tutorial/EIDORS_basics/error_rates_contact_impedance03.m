%Error rate fitting. 
%Assume e = Ch^{p} => log(e) = A + p log(h) and h= C(1/2)^j 
%log(e) = A+log(C) - j(p log(2)) 
%log(err) against jlog(2) will give gradient -p
n_plot=2; %First error 'not asymptotic'so disregard.

%Calculate the best fits for each z_c
for z_cj=1:length(z_c)

for j=1:nmesh
   logh(j) = j*log(2); 
   logcon_func(j)=1; 
   logquad_errL2(j,z_cj) = log(quad_errL2(j,z_cj)); logquad_errH1(j,z_cj) = ...
		 log(quad_errH1(j,z_cj)); logquad_errUM(j,z_cj)=log(quad_errUM(j,z_cj)); 
   logcub_errL2(j,z_cj) = log(cub_errL2(j,z_cj)); logcub_errH1(j,z_cj) = ...
		 log(cub_errH1(j,z_cj)); logcub_errUM(j,z_cj)=log(cub_errUM(j,z_cj));
end

%Extra data point for linear
for j=1:nmesh+1
   loglin_errL2(j,z_cj) = log(lin_errL2(j,z_cj)); loglin_errH1(j,z_cj) = ...
	 		log(lin_errH1(j,z_cj)); loglin_errUM(j,z_cj)=log(lin_errUM(j,z_cj));   
end


%Linear
Alin=[logh',logcon_func']; 
best_linL2(z_cj,:) = (Alin(n_plot:nmesh,:)'*Alin(n_plot:nmesh,:))\...
	Alin(n_plot:nmesh,:)'*loglin_errL2(n_plot:nmesh,z_cj);
best_linH1(z_cj,:) = (Alin(n_plot:nmesh,:)'*Alin(n_plot:nmesh,:))\...
	Alin(n_plot:nmesh,:)'*loglin_errH1(n_plot:nmesh,z_cj);
best_linUM(z_cj,:) = (Alin(n_plot:nmesh,:)'*Alin(n_plot:nmesh,:))\...
	Alin(n_plot:nmesh,:)'*loglin_errUM(n_plot:nmesh,z_cj);

%Quadratic
Aquad=[logh',logcon_func'];
best_quadL2(z_cj,:) = (Aquad(n_plot:nmesh,:)'*Aquad(n_plot:nmesh,:))\...
	Aquad(n_plot:nmesh,:)'*logquad_errL2(n_plot:nmesh,z_cj);
best_quadH1(z_cj,:) = (Aquad(n_plot:nmesh,:)'*Aquad(n_plot:nmesh,:))\...
	Aquad(n_plot:nmesh,:)'*logquad_errH1(n_plot:nmesh,z_cj);
best_quadUM(z_cj,:) = (Aquad(n_plot:nmesh,:)'*Aquad(n_plot:nmesh,:))\...
	Aquad(n_plot:nmesh,:)'*logquad_errUM(n_plot:nmesh,z_cj);

%Cubic
Acub=[logh',logcon_func'];
best_cubL2(z_cj,:) = (Acub(n_plot:nmesh,:)'*Acub(n_plot:nmesh,:))\...
	Acub(n_plot:nmesh,:)'*logcub_errL2(n_plot:nmesh,z_cj);
best_cubH1(z_cj,:) = (Acub(n_plot:nmesh,:)'*Acub(n_plot:nmesh,:))\...
	Acub(n_plot:nmesh,:)'*logcub_errH1(n_plot:nmesh,z_cj);
best_cubUM(z_cj,:) = (Acub(n_plot:nmesh,:)'*Acub(n_plot:nmesh,:))\...
	Acub(n_plot:nmesh,:)'*logcub_errUM(n_plot:nmesh,z_cj);
end

%Plot all linear error rate on one graph against log(z)
figure; hold on;
plot(log10(z_c),-best_linL2(:,1),'ro');
plot(log10(z_c),-best_linH1(:,1),'bx');
plot(log10(z_c),-best_linUM(:,1),'md');
h =legend('L2','H1','UM','Location','NorthWest');
set(h,'FontSize',14);
xlab = xlabel('log(z)'); set(xlab,'FontSize',14);
ylab = ylabel('Rate'); set(ylab,'FontSize',14);      
print_convert error_rates_contact_impedance03a.png

%Plot all quadratic error rate on one graph against log(z)
figure; hold on;
plot(log10(z_c),-best_quadL2(:,1),'ro');
plot(log10(z_c),-best_quadH1(:,1),'bx');
plot(log10(z_c),-best_quadUM(:,1),'md');
h =legend('L2','H1','UM','Location','NorthWest');
set(h,'FontSize',14);
xlab = xlabel('log(z)'); set(xlab,'FontSize',14);
ylab = ylabel('Rate'); set(ylab,'FontSize',14);     
print_convert error_rates_contact_impedance03b.png

%Plot all cubic error rate on one graph against log(z)
figure; hold on;
plot(log10(z_c),-best_cubL2(:,1),'ro');
plot(log10(z_c),-best_cubH1(:,1),'bx');
plot(log10(z_c),-best_cubUM(:,1),'md');
h =legend('L2','H1','I','US','UM','Location','NorthWest');
set(h,'FontSize',14);
xlab = xlabel('log(z)'); set(xlab,'FontSize',14);
ylab = ylabel('Rate'); set(ylab,'FontSize',14);      
print_convert error_rates_contact_impedance03c.png


%Plot all L2 error rates against log(z) for each approx
figure; hold on;
plot(log10(z_c),-best_linL2(:,1),'ro');
plot(log10(z_c),-best_quadL2(:,1),'bx');
plot(log10(z_c),-best_cubL2(:,1),'g+');
h =legend('p=1','p=2','p=3','Location','NorthWest');
set(h,'FontSize',14);
xlab = xlabel('log(z)'); set(xlab,'FontSize',14);
ylab = ylabel('Rate'); set(ylab,'FontSize',14); 
print_convert error_rates_contact_impedance03d.png

%Plot all H1 error rates against log(z) for each approx
figure; hold on;
plot(log10(z_c),-best_linH1(:,1),'ro');
plot(log10(z_c),-best_quadH1(:,1),'bx');
plot(log10(z_c),-best_cubH1(:,1),'g+');
h =legend('p=1','p=2','p=3','Location','NorthWest');
set(h,'FontSize',14);
xlab = xlabel('log(z)'); set(xlab,'FontSize',14);
ylab = ylabel('Rate'); set(ylab,'FontSize',14);      
print_convert error_rates_contact_impedance03e.png

%Plot all UM error rates against log(z) for each approx
figure; hold on;
plot(log10(z_c),-best_linUM(:,1),'ro');
plot(log10(z_c),-best_quadUM(:,1),'bx');
plot(log10(z_c),-best_cubUM(:,1),'g+');
h =legend('p=1','p=2','p=3','Location','NorthWest');
set(h,'FontSize',14);
xlab = xlabel('log(z)'); set(xlab,'FontSize',14);
ylab = ylabel('Rate'); set(ylab,'FontSize',14);      
print_convert error_rates_contact_impedance03f.png
