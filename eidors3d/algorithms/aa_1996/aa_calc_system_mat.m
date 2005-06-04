function s_mat= aa_calc_system_mat( fwd_model, img)
% AA_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.SS  = Unconnected system Matrix
% s_mat.CC  = Connectivity Matrix
% $Id: aa_calc_system_mat.m,v 1.1 2005-06-04 16:39:02 aadler Exp $

p= np_fwd_parameters( fwd_model );

if 0; % old code - less efficient
    SS= spalloc( d*e , d*e, d*d*e);
    for j=1:e
      a=  inv([ ones(d,1) node( :, ELEM(:,j) )' ]);
      SS(d*(j-1)+1:d*j,d*(j-1)+1:d*j)=  ...
           2*a(2:d,:)'*a(2:d,:)/(d-1)/(d-2)/abs(det(a));
    end %for j=1:ELEMs 
end

de= p.d*p.e;
SSiidx= floor([0:de-1]'/p.d)*p.d*ones(1,p.d) + ones(de,1)*(1:p.d) ;
SSjidx= [1:de]'*ones(1,p.d);
SSdata= zeros(de,p.d);
dfact= (p.d-1)*(p.d-2); % Valid for d<=3
for j=1:p.e
  a=  inv([ ones(p.d,1), p.NODE( :, p.ELEM(:,j) )' ]);
  idx= p.d*(j-1)+1:p.d*j;
  SSdata(idx,1:p.d)= 2*a(2:p.d,:)'*a(2:p.d,:)/dfact/abs(det(a));
end %for j=1:ELEMs 
SS= sparse(SSiidx,SSjidx,SSdata);

CC= sparse((1:d*e),ELEM(:),ones(d*e,1), d*e, n);

s_mat.SS = SS;
s_mat.CC = CC;
