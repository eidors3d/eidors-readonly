function data =aa_fwd_solve(R)
% AA_FWD_SOLVE: data= aa_fwd_solve( fwd_model, img)
% Fwd solver for Andy Adler's EIT code
% data = measurements struct
% fwd_model = forward model
% img = image struct

% (C) 1995-2002 Andy Adler
% Ref: Adler & Guardo (1996) IEEE T. Med Imaging
% $Id: aa_fwd_solve.m,v 1.2 2005-06-04 16:39:02 aadler Exp $

p= aa_fwd_parameters( fwd_model );
SS= calc_system_mat( fwd_model, img );

d= size(ELEM,1);        %dimentions+1
e= size(ELEM,2);        %ELEMents
p= size(QQ,2);

if nargin==1
  global NODE SS
  node=NODE;
elseif nargin==2

if 0;
    SS= spalloc( d*e , d*e, d*d*e);
    for j=1:e
      a=  inv([ ones(d,1) node( :, ELEM(:,j) )' ]);
      SS(d*(j-1)+1:d*j,d*(j-1)+1:d*j)=  ...
           2*a(2:d,:)'*a(2:d,:)/(d-1)/(d-2)/abs(det(a));
    end %for j=1:ELEMs 
else
    SSiidx= floor([0:d*e-1]'/d)*d*ones(1,d) + ones(d*e,1)*(1:d) ;
    SSjidx= [1:d*e]'*ones(1,d);
    SSdata= zeros(d*e,d);
    dfact= (d-1)*(d-2); % Note this wont work for d>3
    for j=1:e
      a=  inv([ ones(d,1), node( :, ELEM(:,j) )' ]);
      SSdata(d*(j-1)+1:d*j,1:d)= 2*a(2:d,:)'*a(2:d,:)/dfact/abs(det(a));
    end %for j=1:ELEMs 
    SS= sparse(SSiidx,SSjidx,SSdata);
end

end

if length(R)==max(CPTR)
  R=R(CPTR);
elseif R==0 
  R=zeros(1,size(ELEM,2));
end

i=1:size(SS,2);
d=size(SS,2)/length(R);
z= CC'*SS*sparse(i,i, exp(R(ceil(i/d))) )*CC;

i=2:size(CC,2);

v= [zeros(1,p) ;full( z(i,i)\QQ(i,:) )];
vv= v(MES(1,:),:)- v(MES(2,:),:);
vv= vv(ELS);

