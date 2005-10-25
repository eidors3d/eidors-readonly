function J= aa_e_move_jacobian( fwd_model, img)
% AA_E_MOVE_JACOBIAN: J= aa_e_move_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for EIT, based on conductivity
%   change and movement of electrodes
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc
% $Id: aa_e_move_jacobian.m,v 1.2 2005-10-25 15:49:37 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
stim= fwd_model.stimulation;
s_mat= calc_system_mat( fwd_model, img );

J= conductivity_jacobian( pp, s_mat, stim, fwd_model.gnd_node);
% J= conductivity_jacobian_perturb( pp, 1e-5, img );

d= pp.n_dims+1;
e= pp.n_elem;


% calculate normalized Jacobian
if pp.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,e));
end

function J= conductivity_jacobian( pp, s_mat, stim, gnd_node );
d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
DE= zeros(pp.n_elec, pp.n_stim, e);

idx= 1:pp.n_node;
idx( gnd_node ) = [];
sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat(idx,idx) \ pp.QQ( idx,: );

zinv= zeros(n,n);
zinv(idx,idx) = inv( s_mat(idx,idx) );
zi2E= pp.N2E*zinv;

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
   meas_pat= stim(j).meas_pattern;
   n_meas  = size(meas_pat,1);
   J( idx+(1:n_meas),: ) = meas_pat*DE(:,j,:);
   idx= idx+ n_meas;
end

% Calculate the conductivity jacobian based on a perturbation
% of each element by a delta 
function J= conductivity_jacobian_perturb( pp, delta, img );

J = zeros( pp.n_meas, pp.n_elem );

elem_data = img.elem_data;
d0= fwd_solve( img );
for i=1:pp.n_elem
   img.elem_data   = elem_data;
   img.elem_data(i)= elem_data(i) + delta;
   di= fwd_solve( img );
   J(:,i) = (di.meas - d0.meas) / delta;
end

function J= movement_jacobian( pp, s_mat, stim )
% xy-Movement Jacobian
rr= zeros(1,ne);
for i= 1:2*nm
   idx= elect(ceil(i/2)); % select electrode
   node_xy= rem(i,2) + 1; % select move in x or y
   node= NODE;
   node(node_xy,idx) = node(node_xy,idx) + delta; % electrode movement
   vi = prob_dir( rr, node );
   H(:,ne+i) = (1/delta) * (vi./vh - 1); 
end
