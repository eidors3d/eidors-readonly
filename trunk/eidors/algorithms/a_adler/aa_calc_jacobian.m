function J= aa_calc_jacobian( fwd_model, img)
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

pp= aa_fwd_parameters( fwd_model );
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

FC= aa_system_mat_fields( fwd_model );


if isfield(fwd_model,'coarse2fine')
   DE = jacobian_calc(pp, zi2E, FC, sv, fwd_model.coarse2fine);
   nparam= size(fwd_model.coarse2fine,2);
else
   DE = jacobian_calc(pp, zi2E, FC, sv);
   nparam= e;
end

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
function DE = jacobian_calc(pp, zi2E, FC, sv, c2f);
d= pp.n_dims+1;
dfact= (d-1)*(d-2); % Valid for d<=3

do_c2f = ( nargin==5 );

if ~do_c2f
   DE= zeros(pp.n_elec, pp.n_stim, pp.n_elem);
   zi2E_FCt = zi2E * FC';
   FC_sv   = FC * sv;
   for k= 1:pp.n_elem
       idx= (d-1)*(k-1)+1 : (d-1)*k;
       dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
       DE(:,:,k)= dq;
   end
else
   DE= zeros(pp.n_elec, pp.n_stim, size(c2f,2) );
   if 0 % Code is slower
      zi2E_FCt = zi2E * FC';
      FC_sv   = FC * sv;
      de= pp.n_elem * (d-1);
      for k= 1:size(c2f,2);
          chg_col = kron( c2f(:,k), ones(d-1,1));
          dDD_dEj = spdiags(chg_col,0, de, de);
          dq= zi2E_FCt * dDD_dEj * FC_sv;
          DE(:,:,k)= dq;
      end
   else
      zi2E_FCt = zi2E * FC';
      FC_sv   = FC * sv;
      de= pp.n_elem * (d-1);
      for k= 1:size(c2f,2);
          ff = find( c2f(:,k) );
          lff= length(ff)*(d-1);
          ff1= ones(d-1,1) * ff(:)';
          ffd= (d-1)*ff1 + (-(d-2):0)'*ones(1,length(ff));
          dDD_dEj = spdiags(c2f(ff1,k), 0, lff, lff);
          dq= zi2E_FCt(:,ffd) * dDD_dEj * FC_sv(ffd,:);
          DE(:,:,k)= dq;
      end
   end
end
