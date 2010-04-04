%EidorsDemo1 Demonstrates the use of 2D EIT Package with linear basis
% EidorsDemo1 Demonstrates the use of 2D EIT Package for simulations with linear approksimation basis

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

% Refactored for EIDORS 3.4 by Andy Adler - June 2009
% $Id$


tgt_elems= [374,375,376,601,603,604, ...
            250,254,268,437,449,456];

[fmdl1,fmdl2] = mv_mdl_meshdata;

% Stimulation
stim = mk_stim_patterns(16,1,'{trig}','{ad}',{},1);

% Override with MV's stim pattern
% Trigonometric current pattern.
[jnk,T]=Current(length(fmdl1.electrode),0,'tri');

for i=1:size(T,2)
   stim(i).stim_pattern= T(:,i);
   stim(i).meas_pattern(16,:) = [];
end
%stim(2:end) = [];

fmdl2.stimulation = stim;
fmdl1.stimulation = stim;

Ne2= size(fmdl2.elems,1);
Ne2= size(fmdl2.elems,1);


% Create a sample image
elem_data = 1/400*ones(Ne2,1);
elem_data(tgt_elems) = 1/200;

tgt_img= mk_image(fmdl2,elem_data);

show_fem(tgt_img,[0,1,0])

tgt_img.fwd_model.system_mat = @aa_calc_system_mat;
tgt_img.fwd_model.solve = @aa_fwd_solve;
 meas_aa = fwd_solve( tgt_img );

pp= aa_fwd_parameters( tgt_img.fwd_model );
s_mat= calc_system_mat(tgt_img.fwd_model, tgt_img );
[tgt_img.fwd_model.electrode.z_contact]= deal(50);
v= zeros(pp.n_node,pp.n_stim);

idx= 1:size(s_mat.E,1); idx( tgt_img.fwd_model.gnd_node ) = [];

tol= 1e-5;
v(idx,:)= forward_solver( s_mat.E(idx,idx), pp.QQ(idx,:), tol);


fmdl2.system_mat = @mv_calc_system_mat;
fmdl2.solve = @mv_fwd_solve;
tgt_img.fwd_model= fmdl2;
meas_mv = fwd_solve( tgt_img );

for i=1:10
plot([v(1:1049,i),meas_mv.U.MeasField(:,i)])
pause
end

return

load meshdata % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of element
NNode2=max(size(Node2));                      %The number of nodes
NElement2=max(size(Element2));                %The number of elements

g1=reshape([Node1.Coordinate],2,NNode1)';
H1=reshape([Element1.Topology],3,NElement1)';
g2=reshape([Node2.Coordinate],2,NNode2)';
H2=reshape([Element2.Topology],3,NElement2)';


%disp('Choose a circular inhomogeneity. Left mouse button, center, right button, radius.')
%Ind=ChooseCircle(Node2,Element2);  % Make data for an inhomogeneity.
Ind = tgt_elems;
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=2/400;			  % Conductivity of the inhomogeneity.

clf,Plotinvsol(1./sigma,g2,H2);colorbar,title(['Your resistivity distribution']);drawnow
disp('Press any key to continue...'),pause

disp('Computes the simulated data.')
L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'tri');	  % Trigonometric current pattern.

[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

[U,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,[],'real'); % Simulated data.
Uel=U.Electrode(:);

Agrad1=Agrad*Ind2;   % Group some of the element for the inverse computations


%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.


disp('Solves the full nonlinear inverse problem by regularised Gauss-Newton iteration.')

disp('Initialisations...')

A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);

rho0=Uref.Electrode(:)\U.Electrode(:);

A=UpdateFemMatrix(Agrad,Kb,M,S,1./rho0*ones(size(sigma)));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
Urefel=Uref.Electrode(:);

rho=rho0*ones(size(Agrad1,2),1);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField, ...
           rho,'real');

%Regularisation parameter and matrix

alpha = 0.000005; 
R=MakeRegmatrix(Element1);

iter=5;

disp('Iterations...')

for ii=1:iter
 rho=rho+(J'*J+alpha*R'*R)\(J'*(Uel-Urefel)-alpha*R'*R*rho);
 rhobig=Ind2*rho;
 A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
 Urefel=Uref.Electrode(:);
 J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
 clf,Plotinvsol(rho,g1,H1);colorbar,title([num2str(ii) '. step']);drawnow;
end










