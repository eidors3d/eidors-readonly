%Proposed solutions for the VCIPT Eidors "Hands On" Session
%Banff, Canada:       September 1st 2003 
%Original Solutions:  Nick Polydorides Dept. of Mathematics - UMIST
%Edited by:           David Stephenson Dept. of Chemical Engineering - UMIST


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Introduction:  Setting Up Matlab And Run Demo's.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%These demo loads mesh geometry, attaches electrodes, sets the
%background conductivity, sets the current injection protocol, adds an
%inhomogeneity, calculates the forward and inverse solutions and
%visualises them.

%View the demonstration routine for real conductivity reconstruction by
%copying the black text below into the Matlab command window and pressing
%return.  

demo_real;

%And for complex admittivity reconstruction,

demo_comp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 1:  Loading mesh geometry, generate surface and plot.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the example mesh geometry (vertices matrix and simplices matrix) provided by the datareal MAT file, 
load datareal vtx simp;

%Generate the outer surfaces of the mesh,
[srf] = dubs3(simp);

%Now plot the outer surfaces, 
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3)); 
colormap([0 0 0]); 
daspect([1 1 1]);

%Save workspace as
name_work


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 2:  Set electrode configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the example electrode configuiration provided (then move to "Save")
load datareal elec sels

% |
% |
% |        %Alternatively you may select your own electrode faces using,
% |        [elec_face,sels,cnts,VV] = set_electrodes_n(vtx,srf,[],[],[],[30 60]);
% |        %to select the first electrode face and then,
% |        [elec_face,sels,cnts,VV] = set_electrodes(vtx,srf,elec_face,sels,cnts,VV);
% |        %until you select all the electrode faces.
% |
% |        %Once you have gathered the faces to be assigned as electrodes,  
% |        %provided that you have the same N number of faces per electrode and these are indexed in 
% |        %sequence inside "elec_face", then use something like
% |        elec = reshape(elec_face',(N*3),number of electrodes)';
% |
% |        %i.e for 2 planes of 8 electrodes, of 4 faces each use
% |        elec = reshape(elec_face',(4*3),16)';
% v

%Save workspace as
name_work elec sels

%i.e saving workspace with the new variables "elec" and "sels" defined

%sels is the array holding the indices of the faces selected in 
%elec within the boundary surfaces matrix srf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 3: Select the ground node (or electrode).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%To select the centre node on the top surface of the cylinder as the ground point use,
load datareal gnd_ind
%(This loads the grounding node used by the datareal example)

%Plot the position of the ground node
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3)); 
colormap([0 0 0]); 
daspect([1 1 1]);
hold on;
plot3(vtx(gnd_ind,1),vtx(gnd_ind,2),vtx(gnd_ind,3),'.', ...
                                         'MarkerSize',20,'Color','g');

%Alternatively setting the second electrode elec(2,:) to ground
gnd_ind = elec(2,:) ;

%Save workspace as
name_work gnd_ind 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 4:  Set background conductivity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For a homogeneous isotropic background (real) conductivity of 1 S/m
mat_ref = ones(size(simp,1),1);

%Save workspace as
name_work mat_ref 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 5:  Set current injection pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For some "adjacent" current patterns use the function below 
[I,Ib]=set_3d_currents('{ad}',elec,vtx,gnd_ind,2);

%Alternatively, for some custom made patterns use
Ibt = ones(1,15);
Ibr = - eye(15,15);
Ib = [Ibt;Ibr];
I = [zeros(size(vtx,1),size(Ib,2)); Ib];

%And check that Sum I = 0 (current in = current out) by
sum(I)

%Save workspace as
name_work I Ib 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 6:  Enter contact impedance, compute the system matrix, apply
%reference conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Enter a contact impedance vector, say of 20 Ohms 
zc = 20* ones(size(elec,1),1);

%Compute the E_f matrix, the systrem matrix for 
%the complete electrode model
[Ef] = bld_master_full(vtx,simp,mat_ref,elec,zc);

%Check the size of E_f, and its sparsity
size(Ef)
figure; spy(Ef);

%Then apply the reference conditions
[Er] = ref_master(Ef,gnd_ind);
%for the grounding node gnd_ind or
[Er] = ref_master(Ef,gnd_ind);
%for the grounding electrode gnd_ind

%and then compare numerical ranks before 
rank(full(Ef))
%and after grounding
rank(full(Er))

%Er is also possitive definite, check its eigenvalues/singular values 
eigs(Er)
svds(Ef)

%For the gradient operator D get a percentage measure of its 
%sparsity as
((size(D,1)*size(D,2)) - nnz(D)) / (size(D,1)*size(D,2))*100


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 7:  The reference forward solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Re-compute the system matrix based on mat_ref, call this matrix Eref
[Eref,D,Ela,ppr] = fem_master_full(vtx,simp,mat_ref,gnd_ind,elec,zc,'{n}');

%and calculte the reference forward solution Vref
[Vref] = forward_solver(vtx,Eref,I,1e-5,ppr);

        %or use Vref = Eref\I;

%Then animate the forward solution
potplot(vtx,srf,Vref);

%Save workspace as
name_work zc Eref D Ela Vref 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 8:  Calculate reference measurements, check non-linearity and
%place inhomogeneity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate some reference measurements refH frm Vref
[refH,refV,indH,indV,dfr]=get_3d_meas(elec,vtx,Vref,Ib,2);

%and their index
dfj = dfr(1:2:end);

%Check the non-linearity of the problem (homogeneous isotropic background
%conductivity now set to 3 S/m)
mat_ref3 = 3*ones(size(simp,1),1);

[Eref3,D,Ela,ppr3] = fem_master_full(vtx,simp,mat_ref3,gnd_ind,elec,zc,'{n}');

[Vref3] = forward_solver(vtx,Eref3,I,1e-5,ppr);

[refH3,refV3,indH,indV,dfr3]=get_3d_meas(elec,vtx,Vref3,Ib,2);

%Then plot the difference in the data sets
figure; plot(refH - refH3);

%or alternatively use
figure; plot(refH)
hold on; plot(refH3,'r')

%Create target inhomogeneity of magnitude 1.2 (remember that mat_ref = 1)
[mat,grp] = set_inho(srf,simp,vtx,mat_ref,1.2);

%and maybe a second one as
[mat,grp] = set_inho(srf,simp,vtx,mat,0.85);

%The target conductivity distribution is now mat
%Calculating its corresponding measurements voltageH

[E,D,Ela,pp] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,'{n}');

[V] = forward_solver(vtx,E,I,1e-5,pp);

[voltageH,voltageV,indH,indV,df]=get_3d_meas(elec,vtx,V,Ib,2);

%Plot the perturbation on the data dv due to the introduced inhomogeneities
dv = voltageH - refH;
figure; plot(dv);

%Save workspace as
name_work E D Ela refH voltageH dfj 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 9:  Calculate measurement fields and compute the Jacobian. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The measurement fields v_f for the measurement patterns indH
[v_f] = m_3d_fields(vtx,size(elec,1),indH,Eref,1e-5,gnd_ind);

%Determine the integral of the gradients matrix for first order tetrahedral elements
[IntGrad] = integrofgrad(vtx,simp,mat_ref);

%Define the tolerance
tol = 1e-5

%Compute the Jacobian at the first homogeneous guess mat_ref

[J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,dfj,tol,'{n}');

%Its condition number
cond(J)

%and its singular values
s = svd(J);
figure; semilogy(s,'.');

%Check that
cond(J) >> condest(Eref)

%Are the measurements linearly independent?
rank(J)

%Save workspace as
name_work v_f J 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 10:  Linear inverse solution 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The Tikhonov regularized direct linear solver,

%Add couple of smoothing regularization matrices
%based on neighboring vertices

[R1] = iso_f_smooth(simp,1);
%and faces
[R2] = iso_f_smooth(simp,3);

%Let the relaxation parameter 
a = 1e-7;

%and gengerate some white Gaussian noise
ns = 1e-6*randn(size(refH)); 
%or maybe even more noise

%The corresponding linear regularized solutions for R1 
xt1 = mat_ref + (J.'*J + a*R1.'*R1)\J.' * ((voltageH - refH)+ns);

%and R2
xt2 = mat_ref + (J.'*J + a*R2.'*R2)\J.' * ((voltageH - refH)+ns);

%Write the solutions to a VTK (Visualization ToolKit) file for looking at later,
writevtkfile('xt1',vtx,simp,xt1)

%and
writevtkfile('xt2',vtx,simp,xt2)

%Save workspace as
name_work xt1 xt2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 11:  Linear inverse solution 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The linear Landweber iteration

b = (voltageH - refH)+ns;

x = zeros(size(mat_ref));

t = max(svd(J.'*J));

for i=1:1000 %usually a few thousands are just about enough
   xL = x + t*J.'*(b - J*x);
end

%and the preconditioned conjugate gradient (pcg) solution
[xpcg,flag,relres,iter,resvec] = pcg(J.'*J,J.'*b,1e-4,50);

%Save workspace as
name_work xL xpcg 

%Write the solutions to a VTK (Visualization ToolKit) file using,
writevtkfile('xL',vtx,simp,xL)

%and
writevtkfile('xpcg',vtx,simp,xpcg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 12:  Non-linear inverse solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Newton's nonlinear iterative solver

%For 5 iterations use
[solf] = inv_sol(I,voltageH,1e-5,mat_ref,vtx,simp,elec,2, ...
   zc,'{n}',gnd_ind,a,R1,5);

%Save workspace as
name_work solf

%Write the solutions to a VTK (Visualization ToolKit) file using,
writevtkfile('solf',vtx,simp,solf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 13:  Complex admittivity computations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define a complex reference admittivity
mat_refC = (1+0.5i)*ones(size(simp,1),1);

%and compute the reference system matrix, forward solution, and complex data
[ErefC,D,ElaC,ppr] = fem_master_full(vtx,simp,mat_refC,gnd_ind,elec,zc,'{n}');
[VrefC] = forward_solver(vtx,ErefC,I,1e-5,pp);
[refHC,refVC,indH,indV,df]=get_3d_meas(elec,vtx,VrefC,Ib,2);

%complex measurement fields
[v_fC] = m_3d_fields(vtx,size(elec,1),indH,ErefC,1e-5,gnd_ind);

%Determine the integral of the gradients matrix for first order tetrahedral elements
[IntGradC] = integrofgrad(vtx,simp,mat_refC);

%Define the tolerance
tolC = 1e-5

%complex Jacobian
[JC] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_refC,zc,IntGradC,v_fC,dfj,tolC,'{n}');

%Introduce some complex changes
[matC,grp] = set_inho(srf,simp,vtx,mat_refC,1.1 + 0.6i);
%and 
[matC,grp] = set_inho(srf,simp,vtx,matC,0.9 + 0.6i);

%Compute the data with the inhomogeneity in place
[EC,D,Ela,pp] = fem_master_full(vtx,simp,matC,gnd_ind,elec,zc,'{n}');
[VC] = forward_solver(vtx,EC,I,1e-5,pp);
[voltageHC,voltageVC,indH,indV,df]=get_3d_meas(elec,vtx,VC,Ib,2);

%and a Tikhonov regularized solution
xtC = mat_refC + (JC.'*JC + a*R2.'*R2)\JC.' * ((voltageHC - refHC)+ns);

srC = real(matC); rrC = real(xtC);
siC = imag(matC); riC = imag(xtC);

%Plot the real part of the solution
hsrC = figure;
set(hsrC,'NumberTitle','off');
set(hsrC,'Name','Original conductivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,srC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.63'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,srC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,srC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,srC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,srC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,srC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');

hrrC = figure;
set(hrrC,'NumberTitle','off');
set(hrrC,'Name','Reconstructed conductivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,rrC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.63'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,rrC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,rrC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,rrC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,rrC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,rrC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');

%and the imaginary part

hsiC = figure;
set(hsiC,'NumberTitle','off');
set(hsiC,'Name','Original permitivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,siC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.63'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,siC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,siC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,siC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,siC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,siC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');

hriC = figure;
set(hriC,'NumberTitle','off');
set(hriC,'Name','Reconstructed permitivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,riC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.63'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,riC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,riC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,riC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,riC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,riC,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');

%Save workspace as
name_work complex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 14:  Anisotropic computations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First create an anisotropic admittivity distribution
gamma_xx = 1.0 + 0.5i;
gamma_yy = 1.2 + 0.7i;
gamma_zz = 0.9 + 0.3i;


mat_ani = zeros(size(simp,1)*3,1);
mat_ani(1:3:end) = gamma_xx;
mat_ani(2:3:end) = gamma_yy;
mat_ani(3:3:end) = gamma_zz;

mat_ani = diag(mat_ani);


%Then in bld_master replace the line which calculates Ela 
%....
%This is for the global conductance matrix (includes conductivity)
%Ela = sparse((1:dimen*ns),(1:dimen*ns),kron((a.* Vols).',ones(1,dimen)));
%with

Ela1 = kron( Vols.',ones(1,dimen));

Ela2 = [];
for i=1:ns
Ela2 = sparse(blkdiag(Ela2,a(3*i-2:3*i,3*i-2:3*i)));
end

Ela = diag(Ela1) .* Ela2;

%Save the modified bld_master function under a different name.

%To generate the anisotropic Jacobian the procedure is even more trivial
%All is needed is to comment the line where the x, the y and z components
%of the sensitivity are added together in the isotropic case, i.e.

%Jrow_u = Jrow_x3(1:3:end) + Jrow_x3(2:3:end) + Jrow_x3(3:3:end);

%and use the anisotropic normalized (by the volume) admittivities 
Jrow = Jrow_x3 .* diag(Ela);

%Save workspace as
name_work anisotropic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END.





