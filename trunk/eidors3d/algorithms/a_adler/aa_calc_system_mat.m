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
% $Id: aa_calc_system_mat.m,v 1.13 2007-08-30 03:37:01 aadler Exp $

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
