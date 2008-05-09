function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_calc_system_mat.m,v 1.14 2008-05-09 22:38:00 aadler Exp $

p= aa_fwd_parameters( fwd_model );

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

CC= sparse((1:d*e),p.ELEM(:),ones(d*e,1), d*e, n);

idx= 1:e*d;
elem_sigma = sparse(idx,idx, img.elem_data(ceil(idx/d)) );
s_mat.E= CC'* SS * elem_sigma * CC;


% How to add complete electrode model
bdy = fwd_model.boundary;
for elec= fwd_model.electrode(:)';
   zc=  elec.z_contact;
   bdy_els = zeros(size(bdy,1),1);
   for nd= unique(elec.nodes);
      bdy_els = bdy_els + any(bdy==nd,2);
   end
   ffb = find(bdy_els == size(bdy,2));

   tang_dist=0;
   for ff= ffb(:)'
      bdy_pts= fwd_model.nodes(bdy(ffb,:),:);
      dist=    sqrt( sum(([1,-1]*bdy_pts).^2) );
      tang_dist = tang_dist + dist/zc;

      Ef(m,vr+q) = Ef(m,vr+q) - cali_dist/2 ; % Kv -> Ec  -> Vertical bar
      Ef(n,vr+q) = Ef(n,vr+q) - cali_dist/2 ; % Kv -> Ec
      
      Ef(vr+q,m) = Ef(vr+q,m) - cali_dist/2 ; % Kv' -> Ec' -> Horizontal bar
      Ef(vr+q,n) = Ef(vr+q,n) - cali_dist/2 ; % Kv' -> Ec'
      
      Ef(m,m) = Ef(m,m) + cali_dist/3; % Kz -> E -> Main bar
      Ef(n,n) = Ef(n,n) + cali_dist/3; % Kz -> E
      Ef(m,n) = Ef(m,n) + cali_dist/6; % Kz -> E
      Ef(n,m) = Ef(n,m) + cali_dist/6; % Kz -> E
 

   end
   Ef(vr+q,vr+q) = Ef(vr+q,vr+q) + tang_dist;
   
end

