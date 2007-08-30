function s_mat= dg_calc_system_mat(fwd_model,img)
% DG_CALC_SYSTEM_MAT: SS= aa_calc_system_mat( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: dg_calc_system_mat.m,v 1.3 2007-08-30 03:30:45 aadler Exp $

persistent CC SS;

p= dg_fwd_parameters(fwd_model);
d= p.n_dims+1;
e= p.n_elem;
n= p.n_node;
if isempty(CC) || all(fwd_model.misc.compute_CCandSS=='y')
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
    eidors_msg('dg_calc_system_mat: computing persistent CC and SS',1);
end
idx= 1:e*d;
elem_sigma= sparse(idx,idx, img.elem_data(ceil(idx/d)));
s_mat= CC'*SS*elem_sigma*CC;
