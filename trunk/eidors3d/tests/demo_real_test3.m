function ok= demo_real_test3
% Perform tests based on the demo_real function with new structs

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: demo_real_test3.m,v 1.5 2005-10-30 14:20:30 aadler Exp $

isOctave= exist('OCTAVE_VERSION');

datareal= 'datareal.mat';
datacom=  'datacom.mat';
drt=      'demo_real_test.mat';
if isOctave
    datareal= file_in_loadpath(datareal);
    datacom=  file_in_loadpath(datacom);
    drt    =  file_in_loadpath(drt);
    page_screen_output= 0;
end

% create FEM model structure

load(datareal,'vtx','simp');

demo_mdl= eidors_obj('model', 'Demo real model', ...
                     'nodes', vtx, ...
                     'elems', simp, ...
                     'boundary', dubs3( simp ), ...
                     'solve',      'np_fwd_solve', ...
                     'jacobian',   'np_calc_jacobian', ...
                     'system_mat', 'np_calc_system_mat' );

clear vtx simp

% create FEM model electrodes definitions

load(datareal,'gnd_ind','elec','zc','protocol','no_pl');
perm_sym= '{n}';

demo_mdl= eidors_obj('set', demo_mdl, 'gnd_node', gnd_ind);

for i=1:length(zc)
    demo_mdl.electrode(i).z_contact= zc(i);
    demo_mdl.electrode(i).nodes=     elec(i,:);
end
demo_mdl.misc.perm_sym     = perm_sym;

demo_mdl= eidors_obj('set', demo_mdl);

% create FEM model stimulation and measurement patterns

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

clear gnd_ind elec zc protocol no_pl I Ib
clear indH indV indH_sz meas_pat idx jnk

demo_mdl= eidors_obj('fwd_model', demo_mdl);

% simulate data for homogeneous medium

mat= ones( size(demo_mdl.elems,1) ,1);

homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', mat, ...
                     'fwd_model', demo_mdl );

homg_data=fwd_solve( demo_mdl, homg_img);

% simulate data for inhomogeneous medium

load( datacom ,'A','B') %Indices of the elements to represent the inhomogeneity
mat(A)= mat(A)+0.15;
mat(B)= mat(B)-0.20;

inhomg_img= eidors_obj('image', 'inhomogeneous image', ...
                       'elem_data', mat, ...
                       'fwd_model', demo_mdl );
clear A B mat
inhomg_img = eidors_obj('image', inhomg_img );

inhomg_data=fwd_solve( demo_mdl, inhomg_img);

% create inverse model

% create an inv_model structure of name 'demo_inv'
demo_inv= eidors_obj('inv_model', 'Nick Polydorides EIT inverse', ...
'solve',                  'np_inv_solve', ...
'reconst_type',           'difference', ...
'fwd_model',               demo_mdl);

demo_inv.hyperparameter.value= 1e-8;
demo_inv.image_prior.func= 'np_calc_image_prior';
demo_inv.image_prior.parameters= [3 1];
demo_inv= eidors_obj('set', demo_inv);

% solve inverse model

demo_img= inv_solve( demo_inv, inhomg_data, homg_data);

% verifications

load(drt);

compare_tol( drt.voltageH, inhomg_data.meas, 'voltageH' )
compare_tol( drt.sol, demo_img.elem_data, 'sol' )

J= calc_jacobian( demo_mdl, homg_img );
Jcolsby100=J(:,1:100:size(J,2));
compare_tol( drt.Jcolsby100, Jcolsby100, 'Jcolsby100' )

%Diag_Reg_012= [diag(Reg,0),[diag(Reg,1);0],[diag(Reg,2);0;0]];
%compare_tol( drt.Diag_Reg_012, Diag_Reg_012, 'Diag_Reg_012' )

ok=1;


function compare_tol( cmp1, cmp2, errtext )
% compare matrices and give error if not equal
fprintf(2,'testing parameter: %s ...\n',errtext);

tol= 1e-5;

vd= mean(mean( abs(cmp1 - cmp2) ));
vs= mean(mean( abs(cmp1 + cmp2) ));
if vd/vs > tol
   warning('parameter %s exceeds tolerance %g (=%g)', errtext, tol, vd/vs );
end

