function [org_img, demo_img] = demo_real;
% [org_img, demo_img] = demo_real;
% DEMO to show usage of EIDORS3D
% $Id: demo_real.m,v 1.29 2005-10-12 14:02:38 aadler Exp $
clc;

isOctave= exist('OCTAVE_VERSION');
eidors_msg('log_level',2); % most messages

disp('step 1: create FEM model structure');
%
% create a 'fwd_model' object with name demo_mdl
%

[vtx,simp, bdy] = get_model_elems;

demo_mdl.name = 'demo real model';
demo_mdl.nodes= vtx;
demo_mdl.elems= simp;
demo_mdl.boundary= bdy;
demo_mdl.solve=      'np_fwd_solve';
demo_mdl.jacobian=   'np_calc_jacobian';
demo_mdl.system_mat= 'np_calc_system_mat';

disp('step 2: create FEM model electrodes definitions');

[gnd_ind, electrodes, sym, elec, protocol, no_pl] = get_model_elecs;
demo_mdl.gnd_node=           gnd_ind;
demo_mdl.electrode =         electrodes;
demo_mdl.misc.sym =          sym;

disp('step 3: create FEM model stimulation and measurement patterns');

[stimulations ] = get_model_stim( demo_mdl );
demo_mdl.stimulation= stimulations;

demo_mdl= eidors_obj('fwd_model', demo_mdl); %create object

disp('step 4: simulate data for homogeneous medium');

%
% create a homogeneous image
%

mat= ones( size(demo_mdl.elems,1) ,1);

homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', demo_mdl );

homg_data=fwd_solve( demo_mdl, homg_img);

disp('step 5: simulate data for inhomogeneous medium');
%
% create an inhomogeneous image
% A,B are Indices of the elements to represent the inhomogeneity
%
load( datacom_file ,'A','B')
mat(A)= mat(A)+0.15;
mat(B)= mat(B)-0.20;

inhomg_img= eidors_obj('image', 'inhomogeneous image', ...
                       'elem_data', mat, ...
                       'fwd_model', demo_mdl );

if ~isOctave
%   show_fem( demo_mdl, 1);
    show_fem( inhomg_img , [1 1]);
end


inhomg_data=fwd_solve( demo_mdl, inhomg_img);

disp('step 6: add noise to simulated data');

inhomg_data.meas = inhomg_data.meas + 1e-5*randn(size(inhomg_data.meas));
  homg_data.meas =   homg_data.meas + 1e-5*randn(size(  homg_data.meas));

disp('step 7: create inverse model');

% create an inv_model structure of name 'demo_inv'
demo_inv.name= 'Nick Polydorides EIT inverse';
demo_inv.solve=       'np_inv_solve';
demo_inv.hyperparameter.value = 1e-8;
demo_inv.image_prior.func= 'np_calc_image_prior';
demo_inv.image_prior.parameters= [3 1]; % see iso_f_smooth: deg=1, w=1
demo_inv.reconst_type= 'difference';
demo_inv.fwd_model= demo_mdl;
demo_inv= eidors_obj('inv_model', demo_inv);

disp('step 8: solve inverse model');

demo_img= inv_solve( demo_inv, inhomg_data, homg_data);

disp('step 9: display results');

levels=[ 2.63 2.10 1.72 1.10 0.83 0.10 ];

org_img= eidors_obj('image', 'Simulated inhomogeneities', ...
                    'elem_data', inhomg_img.elem_data - homg_img.elem_data, ...
                    'fwd_model', demo_mdl );
figure; image_levels( org_img, levels );

demo_img.name= 'Reconstructed conductivity distribution';
figure; image_levels( demo_img, levels );

function [vtx,simp,bdy] = get_model_elems;
%bdy : the boundary surfaces (triangles)
%vtx : the vertices of the model (coordinates of the nodes)
%simp: the simplices of the model (connectivity in tetrahedral)
load(datareal_file,'vtx','simp');
bdy= dubs3( simp );



function [gnd_ind, electrodes, sym, elec, protocol, no_pl] = get_model_elecs;
%elec : The electrodes matrix. 
%np_pl : Number of electrode planes (in planar arrangements)
%protocol : Adjacent or Opposite or Customized.
%zc : Contact impedances of the electrodes
%sym : Boolean entry for efficient forward computations 
%sym='{n}';

load(datareal_file,'gnd_ind','elec','zc','protocol','no_pl','sym');

for i=1:length(zc)
    electrodes(i).z_contact= zc(i);
    electrodes(i).nodes=     unique( elec(i,:) );
end

% get the current stimulation patterns
function [stimulations] = get_model_stim( mdl );

load(datareal_file,'protocol','no_pl','sym','elec');
[I,Ib] = set_3d_currents(protocol, ...
                         elec, ...
                         mdl.nodes, ...
                         mdl.gnd_node, ...
                         no_pl);

% get the measurement patterns, only indH is used in this model
%   here we only want to get the meas pattern from 'get_3d_meas',
%   not the voltages, so we enter zeros
[jnk,jnk,indH,indV,jnk] = get_3d_meas( ...
                  elec, mdl.nodes, ...
                  zeros(size(I)), ... % Vfwd
                  Ib, no_pl );
n_elec= size(elec,1);
n_meas= size(indH,1) / size(Ib,2);
for i=1:size(Ib,2)
    stimulations(i).stimulation= 'mA';
    stimulations(i).stim_pattern= Ib(:,i);
    idx= ( 1+ (i-1)*n_meas ):( i*n_meas );
    meas_pat = sparse( (1:n_meas)'*[1,1], ...
                       indH( idx, : ), ...
                       ones(n_meas,2)*[1,0;0,-1], ...
                       n_meas, n_elec );
    stimulations(i).meas_pattern= meas_pat;
end

% Get filename for datareal
function fname= datareal_file;
   fname= 'datareal.mat';
   if exist('OCTAVE_VERSION');
       fname= file_in_loadpath(datareal);
   end

% Get filename for datacom
function fname= datacom_file;
   fname=  'datacom.mat';
   if exist('OCTAVE_VERSION');
       fname=  file_in_loadpath(datacom);
   end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 6.1 R12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
