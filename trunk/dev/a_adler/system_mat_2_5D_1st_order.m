function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: aa_calc_system_mat.m,v 1.5 2006/03/15 21:54:20 aadler Exp $

p= aa_fwd_parameters( fwd_model );

d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

SS = sparse(d*e,d*e);

dfact= factorial(d-1);
for j=1:e
  a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  area = 1/dfact/abs(det(a));
  idx= d*(j-1)+1 : d*j;
  SS(idx,idx)= area*a(2:d,:)'*a(2:d,:);
end %for j=1:ELEMs 

CC= sparse((1:d*e),p.ELEM(:),ones(d*e,1), d*e, n);

idx= 1:e*d;
elem_sigma = sparse(idx,idx, img.elem_data(ceil(idx/d)) );
s_mat.E= CC'* SS * elem_sigma * CC;
