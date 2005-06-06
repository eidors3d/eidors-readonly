function J= aa_calc_jacobian( fwd_model, img)
% AA_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for EIT Alg of Adler & Guardo 1996
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc
% $Id: aa_calc_jacobian.m,v 1.2 2005-06-06 22:35:42 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

% DE_{i,j} is dV_i / dS_j
%  where V_i is change in voltage on electrode i and
%        E_j is change in conductivity on element j
DE= zeros(pp.n_elec, e, pp.n_stim);

idx= 1:pp.n_node;
idx( fwd_model.gnd_node ) = [];
sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat(idx,idx) \ pp.QQ( idx,: );

% connectivity matrix
CC= sparse((1:d*e),pp.ELEM(:),ones(d*e,1), d*e, n);
dfact= (d-1)*(d-2); % Valid for d<=3

for j= 1:e
    a=  inv([ ones(d,1), pp.NODE( :, pp.ELEM(:,j) )' ]);
    dSS_dEj= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));

    idx= d*(j-1)+1 : d*j;
    dq= pp.N2E*(  CC(idx,:)'*dSS_dEj*CC(idx,:)  )*sv;
    DE(:,j,:)= dq;
end

J = zeros( pp.n_meas, e );
idx=0;
for i=1:pp.n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,2);
   J( idx+(1:n_meas),: ) = meas_pat'*DE(:,:,i);
   idx= idx+ n_meas;
end

return;

