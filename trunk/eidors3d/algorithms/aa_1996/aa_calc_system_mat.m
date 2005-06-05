function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix
% $Id: aa_calc_system_mat.m,v 1.3 2005-06-05 13:23:35 aadler Exp $

p= aa_fwd_parameters( fwd_model );

d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

if 0; % old code - less efficient
    SS= spalloc( d*e , d*e, d*d*e);
    for j=1:e
      a=  inv([ ones(d,1) node( :, ELEM(:,j) )' ]);
      SS(d*(j-1)+1:d*j,d*(j-1)+1:d*j)=  ...
           2*a(2:d,:)'*a(2:d,:)/(d-1)/(d-2)/abs(det(a));
    end %for j=1:ELEMs 
end

SSiidx= floor([0:d*e-1]'/d)*d*ones(1,d) + ones(d*e,1)*(1:d) ;
SSjidx= [1:d*e]'*ones(1,d);
SSdata= zeros(d*e,d);
dfact= (d-1)*(d-2); % Valid for d<=3
for j=1:e
  a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  idx= d*(j-1)+1 : d*j;
  SSdata(idx,1:d)= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));
end %for j=1:ELEMs 
SS= sparse(SSiidx,SSjidx,SSdata);

CC= sparse((1:d*e),p.ELEM(:),ones(d*e,1), d*e, n);

idx= 1:e*d;
elem_sigma = sparse(idx,idx, img.elem_data(ceil(idx/d)) );
s_mat= CC'* SS * elem_sigma * CC;
