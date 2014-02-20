function J= test_calc_jacobian( fwd_model, img)
% AA_CALC_JACOBIAN: J= aa_calc_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for current stimulation EIT
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isfield(fwd_model,'coarse2fine'); error('no c2f in this code?'); end

pp= fwd_model_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

idx= 1:size(s_mat.E,1);
idx( fwd_model.gnd_node ) = [];

sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat.E(idx,idx) \ pp.QQ( idx,: );

zi2E= zeros(pp.n_elec, n);
zi2E(:, idx)= pp.N2E(:,idx)/ s_mat.E(idx,idx) ;

FC= system_mat_fields( fwd_model );


DE = jacobian_calc(pp, zi2E, FC, sv);
nparam= e;

J = zeros( pp.n_meas, nparam );
idx=0;
for j= 1:pp.n_stim
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   DEj = reshape( DE(:,j,:), pp.n_elec, nparam );
   J( idx+(1:n_meas),: ) = meas_pat*DEj;
   idx= idx+ n_meas;
end

% calculate normalized Jacobian
if pp.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,nparam));
   
end

% FIXME: The Jacobian calculated is inversed
J= -J;

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
function DE = jacobian_calc(pp, zi2E, FC, sv);
d= pp.n_dims+1;
dfact= (d-1)*(d-2); % Valid for d<=3


zi2E_FCt = zi2E * FC';
FC_sv   = FC * sv;

sz = [pp.n_elec, pp.n_stim, pp.n_elem, d];
if 1
  DE= next_make_c2(sz , zi2E_FCt, FC_sv);
else
  DE= jacobian_loop(sz, zi2E_FCt, FC_sv);
keyboard
end

function DE = next_make_c(sz, zi2E_FCt, FC_sv);
   DE= zeros(sz(1),sz(2),sz(3) );
   d1 = sz(4)-1;
   for k= 1: sz(3);
       idx= d1*(k-1) + (1:d1);
       dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
       DE(:,:,k)= dq;
   end

function DE = next_make_c2(sz, zi2E_FCt, FC_sv);
   DE= zeros(sz(1),sz(2),sz(3) );
   d1 = sz(4)-1;
   FC_sv = FC_sv.';
   for k= 1: sz(3);
       idx= d1*(k-1) + (1:d1);
       M1 =zi2E_FCt(:,idx);
       M2 =   FC_sv(:,idx);
       dq= M1*M2.';
       DE(:,:,k)= dq;
   end
   writefiles(sz, zi2E_FCt, FC_sv, DE)

function writefiles(sz, zi2E_FCt, FC_sv, DE)
save testsave.mat sz zi2E_FCt FC_sv DE

%stupid matlab can't write 3D vectors to binary file
for k = 1:size(DE,3)
   eval(sprintf('DE%02d = DE(:,:,%d);', k, k));
end
clear DE
save testsave.txt -ASCII
