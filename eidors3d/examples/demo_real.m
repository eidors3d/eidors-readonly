% DEMO to show usage of EIDORS3D
% $Id: demo_real.m,v 1.14 2004-07-18 02:46:40 aadler Exp $

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

disp('step 1: create FEM model structure');

load(datareal,'vtx','simp');
%srf : the boundary surfaces (triangles)
%vtx : the vertices of the model (coordinates of the nodes)
%simp: the simplices of the model (connectivity in tetrahedral)

%
% create a 'fwd_model' object with name demo_mdl
%

demo_mdl.type = 'model';
demo_mdl.name = 'demo real model';
demo_mdl.nodes= vtx;
demo_mdl.elems= simp;
demo_mdl.boundary= dubs3( simp );
demo_mdl.solve= 'np_fwd_solve';
demo_mdl.jacobian= 'np_calc_jacobian';

clear vtx simp

disp('step 2: create FEM model electrodes definitions');

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

% TODO: generalize the way that protocol sym no_pl are managed
demo_mdl.misc.sym     = sym;

disp('step 3: create FEM model stimulation and measurement patterns');

% get the current stimulation patterns
[I,Ib] = set_3d_currents(protocol, ...
                         elec, ...
                         demo_mdl.nodes, ...
                         demo_mdl.gnd_node, ...
                         no_pl);
% get the measurement patterns, only indH is used in this model
%   here we only want to get the meas pattern from 'get_3d_meas',
%   not the voltages, so we enter zeros
[jnk,jnk,indH,indV,jnk] = get_3d_meas( ...
                  elec, demo_mdl.nodes, ...
                  zeros(size(I)), ... % Vfwd
                  Ib, no_pl );
n_elec= size(elec,1);
n_meas= size(indH,1) / size(Ib,2);
for i=1:size(Ib,2)
    demo_mdl.stimulation(i).stimulation= 'mA';
    demo_mdl.stimulation(i).stim_pattern= Ib(:,i);
    idx= ( 1+ (i-1)*n_meas ):( i*n_meas );
    meas_pat = sparse( (1:n_meas)'*[1,1], ...
                       indH( idx, : ), ...
                       ones(n_meas,2)*[1,0;0,-1], ...
                       n_meas, n_elec );
    demo_mdl.stimulation(i).meas_pattern= meas_pat;
end

clear gnd_ind elec zc sym protocol no_pl I Ib
clear indH indV indH_sz meas_pat idx jnk


disp('step 4: simulate data for homogeneous medium');

%
% create a homogeneous image
%

homg_img.type = 'image';
homg_img.name = 'homogeneous image';
homg_img.elem_data= ones( size(demo_mdl.elems,1) ,1);
homg_img.fwd_model= demo_mdl;

homg_data=fwd_solve( demo_mdl, homg_img);

disp('step 5: simulate data for inhomogeneous medium');
%
% create an inhomogeneous image
%
mat= ones( size(demo_mdl.elems,1) ,1);
load( datacom ,'A','B') %Indices of the elements to represent the inhomogeneity
mat(A)= mat(A)+0.15;
mat(B)= mat(B)-0.20;

inhomg_img.type = 'image';
inhomg_img.name = 'inhomogeneous image';
inhomg_img.elem_data= mat;
inhomg_img.fwd_model= demo_mdl;
clear A B mat

if ~isOctave
%   show_fem( demo_mdl, 1);
    show_fem( inhomg_img , [1 1]);
end


inhomg_data=fwd_solve( demo_mdl, inhomg_img);

disp('step 6: add noise to simulated data');

% FIXME: I don't understand how this noise should work
% Add noise to data
% dva = homg_data.meas - inhomg_data.meas;
% 
% dc = mean(dva); %DC component of the noise
% noi = dc./7 * ones(length(dva),1) + dc * randn(length(dva),1); %Add the AC component
inhomg_data.meas = inhomg_data.meas + 1e-5*randn(size(inhomg_data.meas));
  homg_data.meas =   homg_data.meas + 1e-5*randn(size(  homg_data.meas));

disp('step 7: create inverse model');

% create an inv_model structure of name 'demo_inv'
demo_inv.name= 'Nick Polydorides EIT inverse';
demo_inv.solve= 'np_inv_solve';
demo_inv.hyperparameter= 1e-8;
demo_inv.type= 'differential';
demo_inv.fwd_model= demo_mdl;

disp('step 8: solve inverse model');

demo_img= inv_solve( demo_inv, inhomg_data, homg_data);

if ~isOctave

    disp('step 9: display results');

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
