function [L2_tot_error,H1semi_tot_error,H1_tot_error,I_err,U_errS,U_errM,U_errSM,timing_solver,DOF]=error_2D_squ_CEM(img,eletype,plot_on)

%Get forward model of the img and the conductivity per element
img.fwd_model.approx_type=eletype;
mdl=img.fwd_model;

%Copy these images for the forward solution
img2=img; mdl2=mdl;


%Modify forward model to ensure elements are reordered and get the new
%nodes and elements
if ~isfield(mdl,'approx_type')    || ...
   strcmp(mdl.approx_type,'tri3') || ...
   strcmp(mdl.approx_type,'tet4')   
else
mdl.approx_type=eletype; [bound,elem,nodes]=fem_1st_to_higher_order(mdl); 
mdl.boundary=bound; mdl.elems=elem; mdl.nodes=nodes;
img.fwd_model=mdl;
end

%Calculate number of nodes and elecs
nnodesF=length(mdl.nodes(:,1)); nelecsF=length(mdl.electrode);
DOF = nnodesF + nelecsF;
tic;
%Calculate stiffness matrix
A_mat = system_mat_higher_order(img);
At=A_mat.E;
ATL = At(1:nnodesF,1:nnodesF);
ATR = At(1:nnodesF,nnodesF+1:nnodesF+nelecsF);
ABL = At(nnodesF+1:nnodesF+nelecsF,1:nnodesF);
ABR = At(nnodesF+1:nnodesF+nelecsF,nnodesF+1:nnodesF+nelecsF);

%Form a newmatirx
AtN=zeros(nnodesF+nelecsF,nnodesF+nelecsF);     
AtN(1:nnodesF,1:nnodesF) = ATL;
AtN(1:nnodesF,nnodesF+1:nnodesF+nelecsF) = 0;
AtN(nnodesF+1:nnodesF+nelecsF,1:nnodesF)=ABL;
AtN(nnodesF+1:nnodesF+nelecsF,nnodesF+1:nnodesF+nelecsF)=-eye(nelecsF);         

%Form a RHS vector
%We do +1 on E1 and -1 on E2
Uvec=[1;-1]; 
RHSvec(1:nnodesF,1) = -ATR*Uvec;
RHSvec(nnodesF+1:nnodesF+nelecsF)=-ABR*Uvec;

%Inverse
uI = AtN \ RHSvec;
volt_nodes=uI(1:nnodesF);
timing_solver=toc;

%Find elemental stiffness matrix on the REFINED MODEL
[~,elemstiff]=mc_calc_stiffness2(mdl,img);


%ANALYTIC Solution

%Okay find the element matrix
elemstruc=mdl.elems; nodestruc=mdl.nodes;

%Find elemental stiffness matrix
[~,elemstiff]=mc_calc_stiffness2(mdl,img);

%Find no. of elements, and initialize vector of H1 errors
nelems=size(elemstruc,1); H1_error=zeros(nelems,1);

%Analytic solution by Fourier decomposition - Two electrodes only
%We have infinite system of matrices U=SA
%U is integral of potential, S is integral of cos products and A are coeffs

%Impedance, COM and half width of each electrode
z_c1=mdl.electrode(1).z_contact; x1=3*pi/10; w1=pi/10;
z_c2=mdl.electrode(2).z_contact; x2=7*pi/10; w2=pi/10; 
%We need ratio of Ul/zl
uz1=Uvec(1)/z_c1; uz2=Uvec(2)/z_c2;

%Number of terms to truncate the analytic and interior terms
n_trunc=1000; n_int_trunc=225;

%
Sc=zeros(n_trunc+1,n_trunc+1);
Uc=zeros(n_trunc+1,1); Ac=zeros(n_trunc+1,1); 

%We first compute the sl(integer) term sl(
sln1=zeros(6*n_trunc+1); sln2=zeros(6*n_trunc+1); 
sln1(3*n_trunc+1)=2*w1; sln2(3*n_trunc+1)=2*w2;
for i=1:3*n_trunc
   sln1(3*n_trunc+1+i) = ( sin(i*(x1+w1))  - sin(i*(x1-w1))  )/i; 
   sln1(3*n_trunc+1-i) = sln1(3*n_trunc+1+i);
   sln2(3*n_trunc+1+i) = ( sin(i*(x2+w2))  - sin(i*(x2-w2))  )/i; 
   sln2(3*n_trunc+1-i) = sln2(3*n_trunc+1+i);   
end

%Assemble matrices and then invert
for i=1:n_trunc+1
   %First workout the t coefficnet
   Uc(i) = uz1*sln1(3*n_trunc+i) + uz2*sln2(3*n_trunc+i) ;   
   for j=1:n_trunc+1       
        Sc(i,j) = 0.5*(1/z_c1)*(sln1(j-i+3*n_trunc+1)+sln1(j+i-2+3*n_trunc+1))...
                + 0.5*(1/z_c2)*(sln2(j-i+3*n_trunc+1)+sln2(j+i-2+3*n_trunc+1));
        if(i==j) %Diagonal term?
            Sc(i,j) = Sc(i,j) + tanh((i-1)*pi)*(i-1)*pi/2;
        end        
   end
end

%Solve and rescle
Ac=Sc\Uc;

%Loop over the nodes and find the boundary
volt_analytic=zeros(length(img.fwd_model.nodes(:,1)),1);
for in=1:length(img.fwd_model.nodes(:,1));
   xin=img.fwd_model.nodes(in,1);
   yin=img.fwd_model.nodes(in,2);
   volt_analytic(in)=Ac(1); %constant
   if(yin==0)
       for k=1:n_trunc
           volt_analytic(in) = volt_analytic(in) + ...
               (Ac(k+1))*cos(k*xin);
       end
   else
       for k=1:n_int_trunc
           volt_analytic(in) = volt_analytic(in) + ...
               (Ac(k+1))*cos(k*xin)*cosh(k*(pi-yin))/cosh(k*pi);
       end
   end
end

%Determine the currents
%Current 1 and 2 constant part
I_analytic(1) = 2*w1/z_c1*(Uvec(1) - Ac(1)); 
I_analytic(2) = 2*w2/z_c2*(Uvec(2) - Ac(1));
for kkk=1:n_trunc
   I_analytic(1) = I_analytic(1) - Ac(kkk+1)/(kkk*z_c1)*( sin(kkk*(x1+w1)) - sin(kkk*(x1-w1)) );
   I_analytic(2) = I_analytic(2) - Ac(kkk+1)/(kkk*z_c2)*( sin(kkk*(x2+w2)) - sin(kkk*(x2-w2)) );   
end
I_FEM(1) = uI(nnodesF+1);
I_FEM(2) = uI(nnodesF+2);

%Calculate the 2-norm error
I_err = norm(I_FEM-I_analytic,2);

%Now solve with curremts
%Ok now lets's solve with the analytic currents
IRHS = zeros(nnodesF+nelecsF,1);
IRHS(nnodesF+1)=I_analytic(1);
IRHS(nnodesF+2)=I_analytic(2);

%Create index vector and eliminate ground node equation from index
groundnode=mdl.gnd_node; idx=1:size(At,1); idx(groundnode)=[];

%Ok now create empty potentials
uU = zeros(nnodesF+nelecsF,1);
uU(idx) = left_divide(At(idx,idx),IRHS(idx,:));
%Subtract zero of potential - x=pi/2 i.e. from symmetry applied potential
uU = uU-0.5*(uU(13)+uU(19));

%Now the potential on stim electrodes
Uvec_FEM(1) = uU(nnodesF+1);
Uvec_FEM(2) = uU(nnodesF+2);

%Error just on stim electrodes
U_errS=norm(Uvec-Uvec_FEM',2);

%Error just on meas (point) electrodes
bound_nodes_not_elecs=img.fwd_model.bound_nodes_not_elecs;
UvecM=volt_analytic(bound_nodes_not_elecs);
UvecM_FEM=uU(bound_nodes_not_elecs);
U_errM = norm(UvecM-UvecM_FEM);

%Error on stim and meas(point) electrodes
UvecSM=[Uvec',UvecM'];
UvecSM_FEM=[Uvec_FEM,UvecM_FEM'];
U_errSM=norm(UvecSM-UvecSM_FEM);

%Plot the solution
if(plot_on==1 && (z_c1==0.00001 || z_c2==1000) && strcmp(mdl.approx_type,'tri10'))
figure; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),volt_nodes,'r*')
hold on; plot3(img.fwd_model.nodes(:,1),img.fwd_model.nodes(:,2),volt_analytic,'b*');
title('Nodal analytic and FEM solutions to CEM');
legend('FEM solution','Analytic solution');
xlabel('x'); ylabel('y'); zlabel('u(x,y)')

if(z_c1==0.00001)
    print_convert error_rates_contact_impedance02a.png;
elseif(z_c1==1000)
    print_convert error_rates_contact_impedance02b.png;
end
end

%Get the elements we want
v=1:nelems;

%Get the basis and gradients at the knot points
eletype=mdl.approx_type; 
if(strcmp(eletype,'tri3'))
    dim=2; order1=0; order2=2;
elseif(strcmp(eletype,'tri6'))    
    dim=2; order1=2; order2=4;
elseif(strcmp(eletype,'tri10'))
    dim=2; order1=4; order2=7;
elseif(strcmp(eletype,'tet4'))
    dim=3; order1=0; order2=2;
elseif(strcmp(eletype,'tet10'))
    dim=3; order1=2; order2=4;
else  
    error('Element type not recognised for integration rules');
end
%Find shape function gradient at knot points
[weight1,xcoord1,ycoord1,zcoord1]=gauss_points(dim,order1);
for kk=1:size(weight1,2)
    dphi_sol(:,:,kk) = element_d_shape_function(eletype,xcoord1(kk),ycoord1(kk),zcoord1(kk));
end
%Find shape functions at knot points (these are higher order ones now)
[weight2,xcoord2,ycoord2,zcoord2]=gauss_points(dim,order2);
for kk=1:size(weight2,2)
    phi_sol(:,kk) = element_shape_function(eletype,xcoord2(kk),ycoord2(kk),zcoord2(kk));
end
%Initialsie the erros
H1semi_error=zeros(nelems,1);
L2_error=zeros(nelems,1);

%Loop through elems and calculate the error
for j=1:length(v)  
    %Get the element
    i=v(j);    
    %Get the node numbers
    eleminodelist=elemstruc(i,:);
    
    %List by row of coordinate on the element
    thise = nodestruc(eleminodelist,:);
    
    %Find the Jacobian of the mapping in 2D and 3D
    jacobianelem = [thise(2,1)-thise(1,1),thise(2,2)-thise(1,2); ...
            thise(3,1)-thise(1,1),thise(3,2)-thise(1,2)];  
    
    %Find the magnitude of the Jacobian of the mapping
    % magjacelem=det(jacobianelem);
    magjacelem=abs(det(jacobianelem));
    
    %L2 error of solution
    %Loop over the weights and evaluate the sensitivty contribution
    for ii=1:size(weight2,2)
        %Find the isoparametric map basis functions at the local point 
        map=element_shape_function('tri3',xcoord2(ii),ycoord2(ii),zcoord2(ii));
        
        %Map our local coordinate to global coordinates using map
        cart(1)=thise(1,1)*map(1)+thise(2,1)*map(2)+thise(3,1)*map(3);
        cart(2)=thise(1,2)*map(1)+thise(2,2)*map(2)+thise(3,2)*map(3);      
        
        %Find the analytic solution at this point  
        analytic_sol = Ac(1);        
        if(cart(2)==0)
           for k=1:n_trunc
           analytic_sol = analytic_sol + (Ac(k+1))*cos(k*cart(1));
           end
        else
           for k=1:n_int_trunc
           analytic_sol = analytic_sol + ...
               (Ac(k+1))*cos(k*cart(1))*cosh(k*(pi-cart(2)))/cosh(k*pi);
           end
        end
        
        %Find the fem solution
        elemi_volt_nodes=volt_nodes(eleminodelist)'; %Vector
        fem_sol=elemi_volt_nodes*phi_sol(:,ii);
        
        %Difference in solution square
        diff_sol=(fem_sol-analytic_sol)^2;
                                
        %Compute the difference and multiply by weight
        diff_sol=diff_sol*weight2(ii);
        
        %Add the contribution to the elemental sensitivity
        L2_error(i)=L2_error(i)+diff_sol;        
    end
    
    
    %H1 error of solution
    %Loop over the weights and evaluate the sensitivty contribution
    for ii=1:size(weight1,2)
        %Find the isoparametric map basis functions at the local point 
        map=element_shape_function('tri3',xcoord1(ii),ycoord1(ii),zcoord1(ii));
        
        %Map our local coordinate to global coordinates using map
        cart(1)=thise(1,1)*map(1)+thise(2,1)*map(2)+thise(3,1)*map(3);
        cart(2)=thise(1,2)*map(1)+thise(2,2)*map(2)+thise(3,2)*map(3);      
                      
        %Find the analytic solution at this point  
        analytic_sol(1) = 0; analytic_sol(2)=0;        
        if(cart(2)==0)
           for k=1:n_trunc
           analytic_sol(1) = analytic_sol(1) - ...
               k*(Ac(k+1))*sin(k*cart(1));
           analytic_sol(2) = analytic_sol(2) - ...
               k*(Ac(k+1))*cos(k*cart(1))*tanh(k*pi);           
           end
        else
           for k=1:n_int_trunc
           analytic_sol(1) = analytic_sol(1) - ...
               k*(Ac(k+1))*sin(k*cart(1))*cosh(k*(pi-cart(2)))/cosh(k*pi);
           analytic_sol(2) = analytic_sol(2) - ...
               k*(Ac(k+1))*cos(k*cart(1))*sinh(k*(pi-cart(2)))/cosh(k*pi);  
           end
        end        
        
        
        %Find the fem solution gradient
        elemi_volt_nodes=volt_nodes(eleminodelist)'; %Vector of local sols
        fem_sol=elemi_volt_nodes*(jacobianelem\dphi_sol(:,:,ii))'; %Gradient
        
        %Difference in solution square
        diff_sol=(fem_sol-analytic_sol)*(fem_sol-analytic_sol)';
                                
        %Compute the difference and multiply by weight
        diff_sol=diff_sol*weight1(ii);
        
        %Add the contribution to the elemental sensitivity
        H1semi_error(i)=H1semi_error(i)+diff_sol;        
    end    
    
    %We have the sensitivity on reference and multiply by Jacobian
    H1semi_error(i)=H1semi_error(i)*magjacelem;   
    L2_error(i)=L2_error(i)*magjacelem; 
            
    %Total error is the sum    
    H1_error(i)=H1semi_error(i)+L2_error(i);
end


%Now find the total H1_error and square root to get norm
H1_tot_error=sum(H1_error); H1_tot_error=sqrt(H1_tot_error);
H1semi_tot_error=sum(H1semi_error); H1semi_tot_error=sqrt(H1semi_tot_error);
L2_tot_error=sum(L2_error); L2_tot_error=sqrt(L2_tot_error);


end