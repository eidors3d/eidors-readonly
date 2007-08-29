function J= aa_calc_jacobian( fwd_model, img)
% AA_CALC_JACOBIAN: J= aa_calc_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for EIT Alg of Adler & Guardo 1996
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_calc_jacobian.m,v 1.1 2007-08-29 09:07:16 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
DE= zeros(pp.n_elec, pp.n_stim, e);

idx= 1:pp.n_node;
idx( fwd_model.gnd_node ) = [];
sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat(idx,idx) \ pp.QQ( idx,: );

   zi2E= zeros(pp.n_elec, n);
%  zi2E(:, idx)= pp.N2E(:,idx)* inv( s_mat(idx,idx) );
   zi2E(:, idx)= pp.N2E(:,idx)/ s_mat(idx,idx) ;

% connectivity matrix
CC= sparse((1:d*e),pp.ELEM(:),ones(d*e,1), d*e, n);
dfact= (d-1)*(d-2); % Valid for d<=3

for k= 1:e
    a=  inv([ ones(d,1), pp.NODE( :, pp.ELEM(:,k) )' ]);
    dSS_dEj= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));

    idx= d*(k-1)+1 : d*k;
    CC_idx = CC(idx,:);
    dq= zi2E * CC_idx' * dSS_dEj * CC_idx * sv;
    DE(:,:,k)= dq;
end

J = zeros( pp.n_meas, e );
idx=0;
for j= 1:pp.n_stim
   meas_pat= fwd_model.stimulation(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   DEj = reshape( DE(:,j,:), pp.n_elec, e );
   J( idx+(1:n_meas),: ) = meas_pat*DEj;
   idx= idx+ n_meas;
end

% calculate normalized Jacobian
if pp.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,e));
   
end

% FIXME: The Jacobian calculated is inversed
J= -J;
