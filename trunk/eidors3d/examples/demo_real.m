% DEMO to show usage of EIDORS3D
% $Id: demo_real.m,v 1.5 2004-07-09 18:51:28 aadler Exp $

clear; 
%clc;

isOctave= exist('OCTAVE_VERSION');

datareal= 'datareal.mat';
datacom=  'datacom.mat';
if isOctave
    datareal= file_in_loadpath(datareal);
    datacom=  file_in_loadpath(datacom);
    page_screen_output= 0;
end

warning('off');
disp('This is a demo for reconstructing conductivity changes')
disp(sprintf('\n'));

load(datareal,'srf','vtx','simp');
%srf : the boundary surfaces (triangles)
%vtx : the vertices of the model (coordinates of the nodes)
%simp: the simplices of the model (connectivity in tetrahedral)

% create a 'fwd_model' object with name demo_mdl

demo_mdl.name = 'demo real model';
demo_mdl.nodes= vtx;
demo_mdl.elems= simp;
demo_mdl.boundary= srf;
demo_mdl.solve= 'np_fwd_solve';

clear srf vtx simp

if ~isOctave
    trimesh(demo_mdl.boundary, ...
            demo_mdl.nodes(:,1), ...
            demo_mdl.nodes(:,2), ...
            demo_mdl.nodes(:,3) );
    axis('image');
    set(gcf,'Colormap',[0 0 0]);
    hold on;
end

disp('This is a cylindrical mesh with homogeneous conductivity distribution of 1')
disp('Wait to attach the electrodes')
disp(sprintf('\n'))

pause(2);

load(datareal,'sels');
%sels :Index in srf matrix denoting the faces to be assigned as electrodes

if ~isOctave
  for u=1:size(sels)
      paint_electrodes(sels(u),srf,vtx);
  end
  hidden('off');
end

clear sels
  
load(datareal,'gnd_ind','elec','zc','protocol','no_pl','sym');
%elec : The electrodes matrix. 
%np_pl : Number of electrode planes (in planar arrangements)
%protocol : Adjacent or Opposite or Customized.
%zc : Contact impedances of the electrodes
%sym : Boolean entry for efficient forward computations 
%sym='{n}';

demo_mdl.gnd_node= gnd_ind;
for i=1:length(zc)
    demo_mdl.electrode(i).z_contact= zc(i);
    demo_mdl.electrode(i).nodes=     elec(i,:);
end

demo_mdl.misc.protocol= protocol;
demo_mdl.misc.sym     = sym;
demo_mdl.misc.no_pl   = no_pl;
% TODO: generalize the way that protocol sym no_pl are managed
clear gnd_ind elec zc sym protocol no_pl

mat_ref = 1*ones(828,1);
data=solve( demo_mdl, mat_ref);
return


elec= zeros(length(demo_mdl.electrode ), ...
            length(demo_mdl.electrode(1).nodes) );
zc  = zeros(length(demo_mdl.electrode ), 1);
for i=1:length(demo_mdl.electrode);
    elec(i,:)= demo_mdl.electrode(i).nodes;
    zc(i)    = demo_mdl.electrode(i).z_contact;
end
[I,Ib] = set_3d_currents(protocol, ...
                         elec, ...
                         demo_mdl.nodes, ...
                         demo_mdl.gnd_node, ...
                         no_pl);

disp('Adjacent current patterns selected') 
disp('Calculating reference measurements')
disp(sprintf('\n'))

pause(2);

mat_ref = 1*ones(828,1); %%%%%%
%Jacobian will be calculated based on this

%Set the tolerance for the forward solver
tol = 1e-5;

[Eref,D,Ela,ppr] = fem_master_full( ...
                demo_mdl.nodes, ...
                demo_mdl.elems, ...
                mat_ref, ...
                demo_mdl.gnd_node, ...
                elec,zc,sym);
[Vref] = forward_solver(demo_mdl.nodes,Eref,I,tol,ppr);
[refH,refV,indH,indV,dfr]=get_3d_meas(elec,demo_mdl.nodes,Vref,Ib,no_pl);
dfr = dfr(1:2:length(dfr)); %Taking just the horrizontal measurements

close;
disp('Allow a local inhomogeneity')
disp(sprintf('\n'))
pause(2);


mat=mat_ref;
load( datacom ,'A','B') %Indices of the elements to represent the inhomogeneity
%figure; [mat,grp] = set_inho(srf,simp,vtx,mat_ref,1.1); 
sA = mat_ref(A(1))+0.15;
sB = mat_ref(B(1))-0.20;
fprintf ('This at %f ',sA)
mat(A) = sA;
mat(B) = sB;

if ~isOctave
figure; 
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3));
axis('image');
set(gcf,'Colormap',[0 0 0]);
hidden('off');
hold on;
repaint_inho(mat,mat_ref,vtx,simp); 
camlight('left');
lighting('flat');
drawnow;

pause(2);
close;
end

disp('Simulating measurements based on ')
disp('the complete electrode model')

[En,D,Ela,ppn] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,sym);
[Vn] = forward_solver(vtx,En,I,tol,ppn,Vref);
[voltageH,voltageV,indH,indV,dfv]=get_3d_meas(elec,vtx,Vn,Ib,no_pl);
dfv = dfv(1:2:end);

if size(dfr)~= size(dfv)
   error('Mismatched measurements')
end


dva = voltageH - refH;
disp('Measurements infused with Gaussian noise ...')


dc = mean(dva); %DC component of the noise
noi = dc./7 * ones(length(dva),1) + dc * randn(length(dva),1); %Add the AC component
dvaG = dva + noi;


disp('Calculating measurement fields')
disp('Please wait, this may take some time ..')
disp(sprintf('\n'))
 [v_f] = m_3d_fields(vtx,32,indH,Eref,tol,gnd_ind);

 disp('Calculating the Jacobian')
[J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,dfr,tol,sym);

disp('Calculating a smoothing prior')
[Reg] = iso_f_smooth(simp,vtx,3,1);

disp('Calculating a linear inverse solution')
disp(sprintf('\n'))

tfac = 1e-8;

sol = (J'*J +  tfac*Reg'*Reg)\J' * dvaG;

if ~isOctave

h1 = figure;
set(h1,'NumberTitle','off');
set(h1,'Name','Simulated inhomogeneities');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,mat-mat_ref,vtx,simp); view(2); grid; colorbar; axis('off'); title('z=2.63'); 
%Calculates also fc. Just once!
subplot(2,3,2); [fc] = slicer_plot_n(2.10,mat-mat_ref,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,mat-mat_ref,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,mat-mat_ref,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,mat-mat_ref,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,mat-mat_ref,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=0.10');


h2 = figure;
set(h2,'NumberTitle','off');
set(h2,'Name','Reconstructed conductivity distribution');
subplot(2,3,1); [fc] = slicer_plot_n(2.63,sol,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=2.63'); 
subplot(2,3,2); [fc] = slicer_plot_n(2.10,sol,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=2.10'); 
subplot(2,3,3); [fc] = slicer_plot_n(1.72,sol,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=1.72'); 
subplot(2,3,4); [fc] = slicer_plot_n(1.10,sol,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=1.10'); 
subplot(2,3,5); [fc] = slicer_plot_n(0.83,sol,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=0.83');
subplot(2,3,6); [fc] = slicer_plot_n(0.10,sol,vtx,simp,fc); view(2); grid; colorbar; axis('off'); title('z=0.10');
disp('Done')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 6.1 R12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
