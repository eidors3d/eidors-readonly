%This demo function shows how the EIT problem can be formulated in a complex
%form. The complex measurements along with a single complex Jacobian are then 
%fed into the reconstruction process. This is a different forulation from the 
%demo_comp.m function.

warning off;
disp('This is a demo for reconstructing admittivity changes of the form a + bi')

load datacom srf vtx simp;
%srf : the boundary surfaces (triangles)
%vtx : the vertices of the model (co-ordinates of the nodes)
%simp: the simplices of the model (connectivity in tetrahedral)

trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
hold on;

disp('This is a cylindrical mesh with homogeneous distribution 1 + 0.5i')
disp('Wait to attach the electrodes')
disp(sprintf('\n'))

pause(2);

load datacom sels;
% FIXME: this needs to use the new show_fem functions
%sels :Index in srf matrix denoting the faces to be assigned as electrodes
% for u=1:size(sels)
%     paint_electrodes(sels(u),srf,vtx);
% end
hidden off;
  
load datacom gnd_ind elec no_pl protocol zc sym;
%gnd_ind : The ground index here a node, can also be an electrode
%elec : The electrodes matrix. 
%np_pl : Number of electrode planes (in planar arrangements)
%protocol : Adjacent or Opposite or Customized.
%zc : Contact impedances of the electrodes
%sym : Boolean entry for efficient forward computations 
%Direct solvers :'{y}' / Iterative : '{n}' 
%sym='{y}';

[I,Ib] = set_3d_currents(protocol,elec,vtx,gnd_ind,no_pl);

disp('Adjacent current patterns selected')
disp('Calculating reference measurements')
disp(sprintf('\n'))

%Set the tolerance for the forward solver
tol = 1e-5;

pause(2);

mat_ref = (1+0.5i)*ones(828,1); %%%%%% 
%Jacobians will be calculated based on this

[Eref,D,Ela,ppr] = fem_master_full(vtx,simp,mat_ref,gnd_ind,elec,zc,sym);
[Vref] = left_divide(Eref,I,tol,ppr);
[refH,refV,indH,indV,dfr]=get_3d_meas(elec,vtx,Vref,Ib,no_pl);
dfr = dfr(1:2:end); %Taking just the horrizontal measurements

close;
disp('Allow a couple of complex changes ...')
disp(sprintf('\n'))
pause(2);


mat=mat_ref;

load datacom A B %Indices of the elements to represent the inhomogeneities
%figure; [mat,grp] = set_inho(srf,simp,vtx,mat_ref,1.2-0.4i); 
sA = 1.2 + 0.4i; % A local complex change or
%sA = 1.2 + 0.5i; %just a local conductivity change
sprintf ('This one  at %f + %f i',real(sA), imag(sA))
mat(A) = sA;

figure; 
subplot(1,2,1);
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
hidden off;
hold on;
repaint_inho(real(mat),real(mat_ref),vtx,simp); title('REAL');
camlight left;
lighting flat;

subplot(1,2,2);
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
hidden off;
hold on;
repaint_inho(imag(mat),imag(mat_ref),vtx,simp); title('IMAGINARY');
camlight left;
lighting flat;

pause(2);
close;
        
%figure; [mat,grp] = set_inho(srf,simp,vtx,mat_ref,0.8-0.6i);       
%sB = 0.8 - 0.6i; % A local complex change
sB = 1 + 0.75i; % or a local permittivity change
sprintf('and this one at  %f + %f i',real(sB), imag(sB))
mat(B) = sB;

figure; subplot(1,2,1); 
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
hidden off;
hold on;
repaint_inho(real(mat),real(mat_ref),vtx,simp); title('REAL');
camlight left;
lighting flat;

subplot(1,2,2); 
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis image;
set(gcf,'Colormap',[0 0 0]);
hidden off;
hold on;
repaint_inho(imag(mat),imag(mat_ref),vtx,simp); title('IMAGINARY');
camlight left;
lighting flat;

pause(2);
close;
        
%[mat,grp] = set_inho(srf,simp,vtx,mat,0.9+0.41i); 
disp(sprintf('\n'))

[En,D,Ela,ppn] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,sym);
[Vn] = left_divide(En,I,tol,ppn);
[voltageH,voltageV,indH,indV,dfv]=get_3d_meas(elec,vtx,Vn,Ib,no_pl);
dfv = dfv(1:2:end);

if size(dfr)~= size(dfv)
   error('Mismatched measurements')
end

disp('Measurements calculated')
disp(sprintf('\n'))

dva = voltageH-refH;

disp('Measurements blended with Gaussian noise ...')
dc = mean(real(dva)); %DC component of the noise
dvrG = dc./7 * ones(length(dva),1) + dc * randn(length(dva),1); %Add the AC component

dc = mean(imag(dva)); %DC component of the noise
dviG = dc./7 * ones(length(dva),1) + dc * randn(length(dva),1); %Add the AC component

dat = (real(dva) + dvrG) + (imag(dva) + dviG)*i;
%dat = dva;

disp('Calculating measurement fields')
[v_f] = m_3d_fields(vtx,size(elec,1),indH,Eref,tol,gnd_ind);

disp('Calculating a single complex jacobian')
[J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,dfv,tol,sym);


disp('Computing a smooting prior')
[Reg] = iso_f_smooth(simp,vtx,1,2);

disp('Calculating a linearised step inverse solution')
disp(sprintf('\n'))
tfac = 1e-7;
sol = (J.'*J +  tfac*Reg'*Reg)\J.' * dat;
sreal = real(sol);
simag = imag(sol);

truereal = real(mat-mat_ref);
trueimag = imag(mat-mat_ref);

v = version;

if str2num(v(1)) > 5

h01 = figure;
set(h01,'NumberTitle','off');
set(h01,'Name','True conductivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,truereal,vtx,simp); view(2); grid; colorbar; axis off; title('z=2.63'); 
%Calculates also fc. Just once!
subplot(2,3,2); [fc] = slicer_plot_n(2.10,truereal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,truereal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,truereal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,truereal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,truereal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');


h02 = figure;
set(h02,'NumberTitle','off');
set(h02,'Name','True scaled permittivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,trueimag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.60'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,trueimag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,trueimag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.70'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,trueimag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,trueimag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.80');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,trueimag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');



h2 = figure;
set(h2,'NumberTitle','off');
set(h2,'Name','Reconstructed conductivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,sreal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.63'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,sreal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,sreal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,sreal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,sreal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,sreal,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');

h3 = figure;
set(h3,'NumberTitle','off');
set(h3,'Name','Reconstructed scaled permittivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,simag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.60'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,simag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,simag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.70'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,simag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,simag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.80');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,simag,vtx,simp,fc); view(2); grid; colorbar; axis off; title('z=0.10');

else
   
   disp('Change the plotting command above to slicer_plot or upgrate to MATLAB 6')
   disp('See demo_real for more details')
   
end

disp('Done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 6.1 R12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
