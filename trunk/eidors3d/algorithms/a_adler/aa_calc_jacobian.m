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
% $Id: aa_calc_jacobian.m,v 1.15 2008-05-09 22:38:00 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

d= pp.n_dims+1;
e= pp.n_elem;
n= pp.n_node;

idx= 1:pp.n_node;
idx( fwd_model.gnd_node ) = [];
sv= zeros(n, pp.n_stim );
sv( idx,:) = s_mat.E(idx,idx) \ pp.QQ( idx,: );

   zi2E= zeros(pp.n_elec, n);
%  zi2E(:, idx)= pp.N2E(:,idx)* inv( s_mat.E(idx,idx) );
   zi2E(:, idx)= pp.N2E(:,idx)/ s_mat.E(idx,idx) ;

% connectivity matrix
CC= sparse((1:d*e),pp.ELEM(:),ones(d*e,1), d*e, n);

if isfield(fwd_model,'coarse2fine')
   DE = jacobian_calc(pp, zi2E, CC, sv, fwd_model.coarse2fine);
   nparam= size(fwd_model.coarse2fine,2);
else
   DE = jacobian_calc(pp, zi2E, CC, sv);
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
   J= J ./ (data.meas(:)*ones(1,e));
   
end

% FIXME: The Jacobian calculated is inversed
J= -J;

% DE_{i,j,k} is dV_i,j / dS_k
%  where V_i is change in voltage on electrode i for
%        stimulation pattern j
%        S_k is change in conductivity on element k
function DE = jacobian_calc(pp, zi2E, CC, sv, c2f);
d= pp.n_dims+1;
dfact= (d-1)*(d-2); % Valid for d<=3

do_c2f = ( nargin==5 );

if do_c2f
   DE= zeros(pp.n_elec, pp.n_stim, size(c2f,2) );
else
   DE= zeros(pp.n_elec, pp.n_stim, pp.n_elem);
end

dSS_dEj= zeros(d,d,pp.n_elem);

for k= 1:pp.n_elem
    a=  inv([ ones(d,1), pp.NODE( :, pp.ELEM(:,k) )' ]);
    F= a(2:d,:) * sqrt(2/dfact/abs(det(a)));
%   dSS_dEj(:,:,k)= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));
    dSS_dEj(:,:,k)= F'*F; % Check transpose for complex?
end

if ~do_c2f
   if 0
      for k= 1:pp.n_elem
          idx= d*(k-1)+1 : d*k;
          CC_idx = CC(idx,:);
          dq= zi2E * CC_idx' * dSS_dEj(:,:,k) * CC_idx * sv;
          DE(:,:,k)= dq;
       end
    elseif 0
      SS= calc_SS(pp);
      FF= calc_FF(pp);
      for k= 1:pp.n_elem
          idx= d*(k-1)+1 : d*k;
          de= pp.n_elem * d;
          dd= zeros(de,1);dd(idx)=1;
          dSS__dEj = SS * spdiags(dd,0, de, de);
          dq= zi2E * CC' * dSS__dEj * CC * sv;
          DE(:,:,k)= dq;
       end
    elseif 0
      FC= calc_FF(pp)*CC;
      zi2E_FCt = zi2E * FC';
      FC_sv   = FC * sv;
      for k= 1:pp.n_elem
          idx= (d-1)*(k-1)+1 : (d-1)*k;
          de= pp.n_elem * (d-1);
          dd= zeros(de,1);dd(idx)=1;
          dDD_dEj = spdiags(dd,0, de, de);
          dq= zi2E_FCt * dDD_dEj * FC_sv;
          DE(:,:,k)= dq;
       end
    else
      FC= calc_FF(pp)*CC;
      zi2E_FCt = zi2E * FC';
      FC_sv   = FC * sv;
      for k= 1:pp.n_elem
          idx= (d-1)*(k-1)+1 : (d-1)*k;
          dq= zi2E_FCt(:,idx) * FC_sv(idx,:);
          DE(:,:,k)= dq;
       end
    end

else
   if 0
      for k= 1:pp.n_elem
          idx= d*(k-1)+1 : d*k;
          CC_idx = CC(idx,:);
          dq= zi2E * CC_idx' * dSS_dEj(:,:,k) * CC_idx * sv;
          chg_col = c2f(k,:);
          for j= find( chg_col );
             DE(:,:,j) = DE(:,:,j) + chg_col(j)*dq;
          end
       end
   elseif 0
      FC= calc_FF(pp)*CC;
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
      FC= calc_FF(pp)*CC;
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

function SS= calc_SS( p) 
d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

SSiidx= floor([0:d*e-1]'/d)*d*ones(1,d) + ones(d*e,1)*(1:d) ;
SSjidx= [1:d*e]'*ones(1,d);
SSdata= zeros(d*e,d);
% FIXME: test_2d_resitor gives wrong results (unless dfact = 4*(d-1))
dfact = (d-2)*(d-1); % to match analytic solution 4*dims
for j=1:e
  a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  idx= d*(j-1)+1 : d*j;
  SSdata(idx,1:d)= 2*a(2:d,:)'*a(2:d,:)/abs(det(a));
end %for j=1:ELEMs 
SS= sparse(SSiidx,SSjidx,SSdata/dfact);

function FF= calc_FF( p) 
d0= p.n_dims+0;
d1= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

FFiidx= floor([0:d0*e-1]'/d0)*d1*ones(1,d1) + ones(d0*e,1)*(1:d1);
FFjidx= [1:d0*e]'*ones(1,d1);
FFdata= zeros(d0*e,d1);
dfact = (d0-1)*d0;
for j=1:e
  a=  inv([ ones(d1,1), p.NODE( :, p.ELEM(:,j) )' ]);
  idx= d0*(j-1)+1 : d0*j;
  FFdata(idx,1:d1)= a(2:d1,:)* sqrt(2/dfact/abs(det(a)));
end %for j=1:ELEMs 
FF= sparse(FFiidx,FFjidx,FFdata)';
