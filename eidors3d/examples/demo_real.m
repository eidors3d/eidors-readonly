% DEMO to show usage of EIDORS3D
% $Id: demo_real.m,v 1.6 2004-07-10 02:40:22 aadler Exp $

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

load(datareal,'sels');
%sels :Index in srf matrix denoting the faces to be assigned as electrodes

if ~isOctave
  for u=1:size(sels)
      paint_electrodes(sels(u),demo_mdl.boundary, ...
                       demo_mdl.nodes);
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

% create a homogeneous image
homg_img.elem_data= ones( size(demo_mdl.elems,1) ,1);
homg_img.fwd_model= demo_mdl;

homg_data=solve( demo_mdl, homg_img);

refH= homg_data.meas;

disp('Allow a local inhomogeneity')
disp(sprintf('\n'))


mat= ones( size(demo_mdl.elems,1) ,1);
load( datacom ,'A','B') %Indices of the elements to represent the inhomogeneity
%figure; [mat,grp] = set_inho(srf,simp,vtx,mat_ref,1.1); 
mat(A)= mat(A)+0.15;
mat(B)= mat(B)+0.15;
clear A B

disp('Simulating measurements based on ')
disp('the complete electrode model')

inhomg_img.elem_data= ones( size(demo_mdl.elems,1) ,1);
inhomg_img.fwd_model= demo_mdl;

if ~isOctave
    figure; 
    trimesh(demo_mdl.boundary, ...
            demo_mdl.nodes(:,1), ...
            demo_mdl.nodes(:,2), ...
            demo_mdl.nodes(:,3) );
    axis('image');
    set(gcf,'Colormap',[0 0 0]);
    hidden('off');
    hold on;
    repaint_inho(inhomg_img.elem_data, ...
                 homg_img.elem_data, ...
                 demo_mdl.nodes, ...
                 demo_mdl.elems); 
    camlight('left');
    lighting('flat');
    drawnow;

    pause(2);
    close;
end

inhomg_data=solve( demo_mdl, inhomg_img);

voltageH= inhomg_data.meas;


dva = voltageH - refH;
disp('Measurements infused with Gaussian noise ...')


% FIXME: this appears to be a strange definition of noise
dc = mean(dva); %DC component of the noise
noi = dc./7 * ones(length(dva),1) + dc * randn(length(dva),1); %Add the AC component
dvaG = dva + noi;

% create an inv_model structure of name 'demo_inv'
demo_inv.name= 'NP EIT inverse';
demo_inv.solve= 'np_inv_solve';
demo_inv.hyperparameter= 1e-8;
demo_inv.type= 'differential';
demo_inv.fwd_model= demo_mdl;

demo_img= inv_solve( demo_inv, homg_data, inhomg_data);

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
