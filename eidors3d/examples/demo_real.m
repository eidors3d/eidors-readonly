% DEMO to show usage of EIDORS3D
% $Id: demo_real.m,v 1.9 2004-07-16 14:36:50 aadler Exp $

clear; 
clc;
warning('off');

isOctave= exist('OCTAVE_VERSION');

datareal= 'datareal.mat';
datacom=  'datacom.mat';
if isOctave
    datareal= file_in_loadpath(datareal);
    datacom=  file_in_loadpath(datacom);
    page_screen_output= 0;
end

disp('step 1: create FEM model');

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

if ~isOctave
   show_fem( demo_mdl)
end
  

% create a homogeneous image
homg_img.elem_data= ones( size(demo_mdl.elems,1) ,1);
homg_img.fwd_model= demo_mdl;

homg_data=fwd_solve( demo_mdl, homg_img);

disp('Allow a local inhomogeneity')
disp(sprintf('\n'))


mat= ones( size(demo_mdl.elems,1) ,1);
load( datacom ,'A','B') %Indices of the elements to represent the inhomogeneity
%figure; [mat,grp] = set_inho(srf,simp,vtx,mat_ref,1.1); 
mat(A)= mat(A)+0.15;
mat(B)= mat(B)-0.20;

inhomg_img.elem_data= mat;
inhomg_img.fwd_model= demo_mdl;
clear A B mat

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

disp('Simulating measurements based on ')
disp('the complete electrode model')

inhomg_data=fwd_solve( demo_mdl, inhomg_img);


% FIXME: I don't understand how this noise should work
% Add noise to data
% dva = homg_data.meas - inhomg_data.meas;
% 
% dc = mean(dva); %DC component of the noise
% noi = dc./7 * ones(length(dva),1) + dc * randn(length(dva),1); %Add the AC component
inhomg_data.meas = inhomg_data.meas + 1e-5*randn(size(inhomg_data.meas));
  homg_data.meas =   homg_data.meas + 1e-5*randn(size(  homg_data.meas));

% create an inv_model structure of name 'demo_inv'
demo_inv.name= 'NP EIT inverse';
demo_inv.solve= 'np_inv_solve';
demo_inv.hyperparameter= 1e-8;
demo_inv.type= 'differential';
demo_inv.fwd_model= demo_mdl;

demo_img= inv_solve( demo_inv, inhomg_data, homg_data);

if ~isOctave

    levels=[ 2.63 2.10 1.72 1.10 0.83 0.10 ];
    org_img= inhomg_img;
    org_img.elem_data= org_img.elem_data - homg_img.elem_data;
    org_img.name= 'Simulated inhomogeneities';
    figure; image_levels( org_img, levels );

    demo_img.name= 'Reconstructed conductivity distribution';
    figure; image_levels( demo_img, levels );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 6.1 R12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
