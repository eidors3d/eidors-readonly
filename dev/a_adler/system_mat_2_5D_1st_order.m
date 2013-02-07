function s_mat= system_mat_2_5D_1st_order( fwd_model, img)
% SYSTEM_MAT_2_5D_1ST_ORDER: 2.5D system matrix
% SS= SYSTEM_MAT_2_5D_1ST_ORDER( FWD_MODEL, IMG)
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix
%
% Parameters
%   fwd_model.system_mat_2_5D_1st_order.k

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: aa_calc_system_mat.m,v 1.5 2006/03/15 21:54:20 aadler Exp $

p= fwd_model_parameters( fwd_model );

try
k= fwd_model.system_mat_2_5D_1st_order.k;
catch
k=0;
end


d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;

SS = sparse(d*e,d*e);

dfact= factorial(d-1);
for j=1:e
  a=  inv([ ones(d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  area = 1/dfact/abs(det(a));
  idx= d*(j-1)+1 : d*j;
  SSe1= a(2:d,:)'*a(2:d,:); 
  SSe2= (ones(3) + eye(3))/12;
  SS(idx,idx)= area*(SSe1 - k^2*SSe2);
end %for j=1:ELEMs 

CC= sparse((1:d*e),p.ELEM(:),ones(d*e,1), d*e, n);

idx= 1:e*d;
elem_sigma = sparse(idx,idx, img.elem_data(ceil(idx/d)) );
s_mat.E= CC'* SS * elem_sigma * CC;
